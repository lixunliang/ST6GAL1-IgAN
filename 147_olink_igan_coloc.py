#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
147_olink_igan_coloc.py  —  IgA肾病候选蛋白：Olink全区域共定位（2026-10-07）

145步只用了每个蛋白的领头pQTL。这里用146步下载的UKB-PPP全量数据，对25个候选蛋白做：
  1 领头cis-pQTL上下游500 kb内，蛋白信号 × IgA肾病信号的共定位（PP3=两个不同信号，PP4=同一个信号）
    分别对欧洲、亚洲、合并、去Immunochip四份IgA肾病数据
  2 两边的最强位点互查：蛋白最强位点在IgA肾病里的P；IgA肾病最强位点在蛋白里的z；两者相距多远
  3 同一位点上相邻基因的蛋白互相对照（例如16号染色体上的ITGAM与PYDC1）
位点按 染色体:GRCh37坐标 + 等位基因 匹配（Olink的ID列自带GRCh37坐标）。
注意：pQTL来自欧洲人群，与亚洲IgA肾病做共定位时两边LD结构不同，PP4会被低估，亚洲一列只作参考。
先验与前面各步相同（p1=p2=1e-4，p12=1e-5，W=0.15²）。

不联网。约20分钟（读25个tar的各一条染色体 + 扫4个IgA肾病文件）。
输出  results\\olink_igan_coloc.tsv ; logs\\olink_igan_coloc_report.txt
"""
import csv
import glob
import gzip
import importlib.util
import io
import math
import os
import re
import sys
import tarfile
import time
import traceback

import numpy as np

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "olink_igan_coloc_report.txt")
INSTR = os.path.join(EXT, "olink_cis_instruments.tsv")
TAR_DIR = os.path.join(PROJ, "raw", "protein", "ukbppp")
K = os.path.join(PROJ, "raw", "outcome", "igan_kiryluk2023")
IGAN = [("欧洲", os.path.join(K, "IgAN_Combined_metaanalysis_European_only.txt")),
        ("亚洲", os.path.join(K, "IgAN_Combined_metaanalysis_Asian_only.txt")),
        ("合并", os.path.join(K, "IgAN_Combined_metaanalysis.txt")),
        ("去IC", os.path.join(K, "IgAN_Combined_metaanalysis_without_immunochip.txt"))]
GENES = ["TNFSF13", "TNFRSF13B", "TNFRSF13C", "CFH", "CFHR5", "CFB", "C3", "ST6GAL1", "ITGAM", "PYDC1", "FCAR",
         "FCGR2A", "FCGR2B", "FCRL2", "FCRL3", "PADI4", "EIF2AK3", "NFKB1", "IL7R", "CD28", "RNASET2", "SAT2",
         "TNFRSF6B", "CXCL6", "PIGR"]
if os.environ.get("MR_TEST_GENES"):
    GENES = os.environ["MR_TEST_GENES"].split(",")
# 领头位点的注释（来自Sun 2023 附表ST9），提示可能的表位效应或交叉反应
NOTE = {"TNFSF13": "领头位点为TNFSF13错义变异", "TNFRSF13B": "领头位点为TNFRSF13B错义变异",
        "FCGR2B": "领头位点是FCGR2A的错义变异rs1801274（H131R），可能是抗体交叉反应", "FCRL3": "领头位点为FCRL3错义变异",
        "IL7R": "领头位点为IL7R错义变异", "SAT2": "领头位点为SAT2错义变异", "CFB": "领头位点为CFB错义变异，且在MHC区",
        "PIGR": "领头位点为PIGR错义变异（频率2.4%）", "CFHR5": "领头位点远离CFHR5基因（在KCNT2内，频率1.7%）"}
WIN = 500_000
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
W, P1, P2, P12 = 0.15 ** 2, 1e-4, 1e-4, 1e-5
_out, _rows = [], []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def coloc(b1, s1, b2, s2):
    def labf(b, s):
        r = W / (W + s ** 2)
        return 0.5 * (np.log(1 - r) + r * (b / s) ** 2)

    def lse(x):
        m = np.max(x)
        return m + np.log(np.sum(np.exp(x - m)))
    l1, l2 = labf(b1, s1), labf(b2, s2)
    L1, L2, L4 = lse(l1), lse(l2), lse(l1 + l2)
    L3 = L1 + L2 + math.log1p(-math.exp(min(L4 - L1 - L2, -1e-12)))
    lh = np.array([0.0, math.log(P1) + L1, math.log(P2) + L2, math.log(P1) + math.log(P2) + L3, math.log(P12) + L4])
    return np.exp(lh - lse(lh))


def read_olink(tar_path, chrom, lo, hi):
    """返回 {pos37: (a0, a1, beta, se)}（a1为效应等位基因；只留单碱基位点；同坐标多等位时留|z|大的）"""
    out = {}
    with tarfile.open(tar_path, "r") as tf:
        cand = [m for m in tf.getmembers() if m.isfile() and re.search(rf"chr{chrom}[^0-9]", m.name, re.I)]
        if not cand:
            return out
        fobj = tf.extractfile(cand[0])
        raw = gzip.GzipFile(fileobj=fobj) if cand[0].name.endswith(".gz") else fobj
        st = io.TextIOWrapper(io.BufferedReader(raw, 1 << 20), encoding="utf-8", errors="replace")
        ix = {h.upper(): i for i, h in enumerate(st.readline().split())}
        for line in st:
            v = line.split()
            try:
                q = v[ix["ID"]].split(":")
                pos = int(q[1])
                if pos < lo or pos > hi:
                    continue
                a0, a1 = v[ix["ALLELE0"]].upper(), v[ix["ALLELE1"]].upper()
                if len(a0) != 1 or len(a1) != 1:
                    continue
                b, se = float(v[ix["BETA"]]), float(v[ix["SE"]])
            except (ValueError, IndexError, KeyError):
                continue
            if se > 0 and (pos not in out or abs(b / se) > abs(out[pos][2] / out[pos][3])):
                out[pos] = (a0, a1, b, se)
    return out


def scan_igan(path, regions):
    """regions: {chr: [(lo, hi)]}；返回 {(chr,pos): (a1, a2, beta, se)}"""
    d, n = {}, 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        low = [h.lower() for h in f.readline().split()]
        ic, ip, ia1, ia2, ib, ise = (low.index(x) for x in ("chr", "bp_hg19", "a1", "a2", "beta", "se"))
        for line in f:
            n += 1
            v = line.split()
            iv = regions.get(v[ic])
            if iv is None:
                continue
            try:
                pos = int(v[ip])
            except ValueError:
                continue
            if not any(lo <= pos <= hi for lo, hi in iv):
                continue
            try:
                b, se = float(v[ib]), float(v[ise])
            except ValueError:
                continue
            if se > 0 and math.isfinite(b):
                d[(v[ic], pos)] = (v[ia1].upper(), v[ia2].upper(), b, se)
    return d, n


def sgn(a1, a0, ea, oa):
    """Olink效应等位基因a1/a0 对齐到 IgAN的ea/oa"""
    if {a1, a0} in ({"A", "T"}, {"C", "G"}):
        return 0
    if (a1, a0) == (ea, oa):
        return 1
    if (a1, a0) == (oa, ea):
        return -1
    c = (COMP.get(ea, "?"), COMP.get(oa, "?"))
    if (a1, a0) == c:
        return 1
    if (a1, a0) == (c[1], c[0]):
        return -1
    return 0


def run():
    lead = {}
    with open(INSTR, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["source"] == "ST9_lead" and r["Gene"] in GENES and r["Gene"] not in lead:
                lead[r["Gene"]] = (r["CHR"], int(r["POS37"]), r["rsid"])
    olink, regions = {}, {}
    log("A. 读Olink全量数据（领头cis-pQTL ±500 kb）")
    for g in GENES:
        if g not in lead:
            log(f"  {g}：工具变量表里没有领头位点，跳过")
            continue
        tars = [t for t in glob.glob(os.path.join(TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith(g + "_")]
        if not tars:
            log(f"  {g}：没有tar文件，跳过")
            continue
        c, p, rs = lead[g]
        t0 = time.time()
        try:
            d = read_olink(tars[0], c, p - WIN, p + WIN)
        except Exception as e:  # noqa: BLE001
            log(f"  {g}：读tar失败 {type(e).__name__}: {e}")
            continue
        if not d:
            log(f"  {g}：窗口内没读到位点")
            continue
        olink[g] = d
        regions.setdefault(c, []).append((p - WIN, p + WIN))
        top = max(d.items(), key=lambda kv: abs(kv[1][2] / kv[1][3]))
        log(f"  {g:<10} chr{c}:{p - WIN:,}-{p + WIN:,}  位点 {len(d):,}  最强pQTL chr{c}:{top[0]:,} z={top[1][2] / top[1][3]:+.1f}"
            f"  （领头 {rs}）  {time.time() - t0:.0f}s")
    log("")
    log("B. 扫IgA肾病文件")
    ig = {}
    for name, path in IGAN:
        if not os.path.exists(path):
            log(f"  {name}：文件不存在")
            continue
        d, n = scan_igan(path, regions)
        ig[name] = d
        log(f"  {name}：扫描 {n:,} 行，窗口内 {len(d):,} 个位点")
    log("")
    log("C. 共定位。每个蛋白一段；'互查'两行说明两边最强位点是否是同一处")
    summ = {}
    for g in GENES:
        if g not in olink:
            continue
        c, p, rs = lead[g]
        d = olink[g]
        log(f"  {g}" + (f"  [{NOTE[g]}]" if g in NOTE else ""))
        for name, _p in IGAN:
            idat = ig.get(name)
            if not idat:
                continue
            pos, b1, s1, b2, s2 = [], [], [], [], []
            for q, (a0, a1, b, se) in d.items():
                o = idat.get((c, q))
                if o is None:
                    continue
                s = sgn(a1, a0, o[0], o[1])
                if s == 0:
                    continue
                pos.append(q); b1.append(b); s1.append(se); b2.append(s * o[2]); s2.append(o[3])
            if len(pos) < 50:
                log(f"    {name}：共有位点只有 {len(pos)} 个，跳过")
                continue
            pos = np.array(pos); b1 = np.array(b1); s1 = np.array(s1); b2 = np.array(b2); s2 = np.array(s2)
            pp = coloc(b1, s1, b2, s2)
            z1, z2 = b1 / s1, b2 / s2
            i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
            wald = b2[i1] / b1[i1]
            wse = abs(s2[i1] / b1[i1])
            log(f"    {name}  位点 {len(pos):>5}  PP3={pp[3]:.2f} PP4={pp[4]:.2f}  "
                f"蛋白最强位点 {c}:{pos[i1]:,}（z={z1[i1]:+.0f}）在IgA肾病 P={pz(z2[i1]):.1e}，Wald={wald:+.2f}；"
                f"IgA肾病最强位点 {c}:{pos[i2]:,}（P={pz(z2[i2]):.1e}）在蛋白 z={z1[i2]:+.1f}；相距 {abs(int(pos[i1]) - int(pos[i2])) / 1000:.0f} kb")
            _rows.append([g, name, len(pos), f"{pp[3]:.4g}", f"{pp[4]:.4g}", f"{c}:{pos[i1]}", f"{z1[i1]:.2f}",
                          f"{pz(z2[i1]):.3g}", f"{wald:.4g}", f"{wse:.4g}", f"{c}:{pos[i2]}", f"{pz(z2[i2]):.3g}",
                          f"{z1[i2]:.2f}", abs(int(pos[i1]) - int(pos[i2])), NOTE.get(g, "")])
            summ[(g, name)] = (pp[4], pp[3], wald, pz(z2[i1]), pz(z2[i2]))
    log("")
    log("D. 汇总（PP4：欧洲 / 亚洲 / 合并 / 去IC；后面是合并数据里IgA肾病在该窗口的最小P）")
    for g in GENES:
        if not any((g, n) in summ for n, _p in IGAN):
            continue
        cells = [f"{summ[(g, n)][0]:.2f}" if (g, n) in summ else " -  " for n, _p in IGAN]
        mp = summ.get((g, "合并"), (0, 0, 0, 1, float("nan")))[4]
        verdict = "共定位成立" if max((summ[(g, n)][0] for n, _p in IGAN if (g, n) in summ), default=0) >= 0.8 else \
            ("两个不同信号" if max((summ[(g, n)][1] for n, _p in IGAN if (g, n) in summ), default=0) >= 0.8 else "不确定/IgA肾病信号弱")
        log(f"    {g:<10} " + " / ".join(cells) + f"   窗口内IgA肾病最小P={mp:.1e}   {verdict}")
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "olink_igan_coloc.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("Gene\tIgAN_set\tn_snps\tPP3\tPP4\ttop_pQTL\ttop_pQTL_z\ttop_pQTL_IgAN_p\twald\twald_se\ttop_IgAN\t"
                "top_IgAN_p\ttop_IgAN_pQTL_z\tdistance_bp\tnote\n")
        for r in _rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(LOG_DIR, exist_ok=True)
    try:
        run()
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")


if __name__ == "__main__":
    main()
