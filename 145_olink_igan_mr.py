#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
145_olink_igan_mr.py  —  Olink cis-pQTL → IgA nephropathy screen

1 Olink（UKB-PPP，Sun 2023附表ST9/ST16）的cis工具变量 → IgA肾病（Kiryluk 2023 欧洲/亚洲/合并/去Immunochip）
  及FinnGen 5个肾脏病种。工具变量表 extracted\\olink_cis_instruments.tsv 由附表整理而来（约1,930个蛋白，只保留单碱基位点）：
  每个蛋白优先用ST9的领头cis-pQTL；领头位点在结局里找不到时，退到ST16里最强的独立cis信号
  （ST16的效应是标准化基因型尺度，已换算成每个等位基因）。单位点Wald比值。
  Kiryluk文件按 染色体:GRCh37坐标 匹配，FinnGen按 染色体:GRCh38坐标 匹配，并核对等位基因。
  注意：只有领头位点、没有全区域数据，所以这一步没有共定位，只是筛选。
不联网。约15分钟（Kiryluk 4个文件各扫一遍，FinnGen 5个文件各扫一遍）。
输出  results\\olink_kidney_disease_mr.tsv ; logs\\olink_igan_mr_report.txt
"""
import csv
import gzip
import importlib.util
import math
import os
import sys
import time
import traceback
from collections import Counter

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "olink_igan_mr_report.txt")
INSTR = os.path.join(EXT, "olink_cis_instruments.tsv")
CACHE = os.path.join(EXT, "kd_olink_instrument_assoc.tsv")
SOMA = os.path.join(RES, "kidney_disease_protein_mr.tsv")
SOMA_COLOC = os.path.join(RES, "kidney_disease_coloc.tsv")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
IG = ["IGAN_EUR", "IGAN_ASN", "IGAN_ALL", "IGAN_NOIC", "IGAN_FG"]
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def align(ea, oa, ea2, oa2):
    if (ea, oa) == (ea2, oa2):
        return 1
    if (ea, oa) == (oa2, ea2):
        return -1
    c = (COMP.get(ea2, "?"), COMP.get(oa2, "?"))
    if (ea, oa) == c:
        return 1
    if (ea, oa) == (c[1], c[0]):
        return -1
    return 0


def scan(m144, path, fmt, keys):
    """keys: set of (chr, pos)。返回 {(chr,pos): [(ea, oa, beta, se), ...]}"""
    d, n = {}, 0
    with m144.opener(path) as f:
        low = [h.lower().lstrip("#") for h in f.readline().split()]
        if fmt == "kiryluk":
            ic, ip, ia1, ia2, ib, ise = (low.index(x) for x in ("chr", "bp_hg19", "a1", "a2", "beta", "se"))
        else:
            ic, ip, ia1, ia2, ib, ise = (low.index(x) for x in ("chrom", "pos", "alt", "ref", "beta", "sebeta"))
        mx = max(ic, ip, ia1, ia2, ib, ise)
        for line in f:
            n += 1
            v = line.split("\t") if "\t" in line else line.split()
            if len(v) <= mx:
                continue
            k = (v[ic].replace("chr", ""), v[ip])
            if k not in keys:
                continue
            try:
                b, se = float(v[ib]), float(v[ise])
            except ValueError:
                continue
            if se > 0 and math.isfinite(b):
                d.setdefault(k, []).append((v[ia1].strip().upper(), v[ia2].strip().upper(), b, se))
    return d, n


def bh(ps):
    n = len(ps)
    o = sorted(range(n), key=lambda i: ps[i])
    out, m = [0.0] * n, 1.0
    for rank in range(n, 0, -1):
        i = o[rank - 1]
        m = min(m, ps[i] * n / rank)
        out[i] = m
    return out


def run():
    spec = importlib.util.spec_from_file_location("_m144", os.path.join(HERE, "144_proteome_to_kidney_diseases.py"))
    m144 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m144)
    # ---------- 1 Olink工具变量
    prot = {}
    with open(INSTR, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b, se = fnum(r["Beta"]), fnum(r["SE"])
            if not (se > 0 and b == b and b != 0):
                continue
            prot.setdefault(r["ProteinID"], []).append(dict(
                gene=r["Gene"], src=r["source"], rs=r["rsid"], c=r["CHR"], p37=r["POS37"], p38=r["POS38"],
                ea=r["EA"].upper(), oa=r["OA"].upper(), b=b, se=se, lp=fnum(r["log10p"]), mhc=int(fnum(r["MHC"]) or 0)))
    for v in prot.values():
        v.sort(key=lambda x: (x["src"] != "ST9_lead", -x["lp"]))
    log(f"Olink蛋白 {len(prot):,} 个，cis工具变量 {sum(len(v) for v in prot.values()):,} 个")
    keys37 = {(x["c"], x["p37"]) for v in prot.values() for x in v}
    keys38 = {(x["c"], x["p38"]) for v in prot.values() for x in v if x["p38"]}
    assoc = {}
    if os.path.exists(CACHE):
        with open(CACHE, "r", encoding="utf-8") as f:
            f.readline()
            for line in f:
                o, c, p, ea, oa, b, se = line.rstrip("\n").split("\t")
                assoc.setdefault(o, {}).setdefault((c, p), []).append((ea, oa, float(b), float(se)))
    log("")
    log("A. 抽工具变量")
    for o, lab, path, fmt in m144.OUTFILES:
        if o in assoc:
            log(f"  {o:<10} 已有缓存 {len(assoc[o]):,} 个位点")
            continue
        if not os.path.exists(path):
            log(f"  {o:<10} !! 文件不存在：{path}")
            continue
        t0 = time.time()
        d, n = scan(m144, path, fmt, keys37 if fmt == "kiryluk" else keys38)
        assoc[o] = d
        log(f"  {o:<10} 扫描 {n:,} 行，抽到 {len(d):,} 个位点（{(time.time() - t0) / 60:.1f} 分钟）")
        with open(CACHE, "w", encoding="utf-8", newline="") as f:
            f.write("outcome\tchr\tpos\tEA\tOA\tbeta\tse\n")
            for oo, dd in assoc.items():
                for (c, p), lst in dd.items():
                    for ea, oa, b, se in lst:
                        f.write(f"{oo}\t{c}\t{p}\t{ea}\t{oa}\t{b:.8g}\t{se:.8g}\n")
    # ---------- 2 Wald
    log("")
    log("B. Olink cis-MR（Wald；每SD蛋白的ln OR）")
    res, rows_all = {}, []
    for o, lab, _p, fmt in m144.OUTFILES:
        d = assoc.get(o)
        if not d:
            continue
        rows = []
        for pid, ins in prot.items():
            for x in ins:
                k = (x["c"], x["p37"]) if fmt == "kiryluk" else (x["c"], x["p38"])
                hit = None
                for ea, oa, b, se in d.get(k, []):
                    s = align(x["ea"], x["oa"], ea, oa)
                    if s and {x["ea"], x["oa"]} not in ({"A", "T"}, {"C", "G"}):
                        hit = (s * b, se)
                        break
                if hit is None:
                    continue
                w, wse = hit[0] / x["b"], abs(hit[1] / x["b"])
                rows.append(dict(outcome=o, ProteinID=pid, Gene=x["gene"], source=x["src"], rsid=x["rs"],
                                 chr=x["c"], pos37=x["p37"], pqtl_beta=x["b"], pqtl_log10p=x["lp"],
                                 beta=w, se=wse, p=pz(w / wse), MHC=x["mhc"]))
                break
        if not rows:
            continue
        q = bh([r["p"] for r in rows])
        for r, qv in zip(rows, q):
            r["FDR"] = qv
            r["bonferroni"] = int(r["p"] < 0.05 / len(rows))
            res[(o, r["Gene"])] = r
        rows.sort(key=lambda r: r["p"])
        rows_all += rows
        nb = [r for r in rows if r["bonferroni"]]
        log(f"  {o:<10} {lab}：可检验 {len(rows):,}（用领头位点 {sum(r['source'] == 'ST9_lead' for r in rows):,}）；"
            f"Bonferroni {len(nb)}（MHC区 {sum(r['MHC'] for r in nb)}）；FDR<0.05 {sum(r['FDR'] < 0.05 for r in rows)}")
        show = [r for r in rows if not r["MHC"]][:15] if o.startswith("IGAN") else rows[:8]
        for r in show:
            log(f"      {r['Gene']:<10} OR={math.exp(max(min(r['beta'], 20), -20)):.2f} beta={r['beta']:+.3f} P={r['p']:.1e} "
                f"FDR={r['FDR']:.1e} {r['rsid']} pQTL={r['pqtl_beta']:+.2f}{'' if r['source'] == 'ST9_lead' else ' (次级信号)'}"
                f"{'  [MHC]' if r['MHC'] else ''}")
    cols = ["outcome", "ProteinID", "Gene", "source", "rsid", "chr", "pos37", "pqtl_beta", "pqtl_log10p", "beta", "se",
            "p", "FDR", "bonferroni", "MHC"]
    with open(os.path.join(RES, "olink_kidney_disease_mr.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows_all:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    log("")
    log("C. IgA肾病通路的重点蛋白（Olink）：欧洲 / 亚洲 / 合并")
    for gname in ["TNFSF13", "TNFSF13B", "TNFRSF13B", "TNFRSF17", "FCAR", "PIGR", "C3", "CFHR5", "CFH", "CFB", "MASP1",
                  "ITGAM", "OSM", "IL6R", "ST6GAL1", "TGFB1", "FCRL3", "VAV3", "CD79B", "TNFSF8", "CCL20", "IL2RA",
                  "CD28", "GZMB", "PSAP", "AGT", "ACE"]:
        cells = []
        for o in ("IGAN_EUR", "IGAN_ASN", "IGAN_ALL"):
            r = res.get((o, gname))
            cells.append(f"{r['beta']:+.2f}({r['p']:.0e})".ljust(15) if r else "-".ljust(15))
        log(f"    {gname:<10} " + " ".join(cells))


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    for d in (LOG_DIR, RES):
        os.makedirs(d, exist_ok=True)
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
