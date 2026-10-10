#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
144_proteome_to_kidney_diseases.py  —  循环蛋白 → 具体肾脏疾病（首选IgA肾病）（2026-10-07）

暴露  D:\\肥胖国自然\\extracted\\cis_instruments.tsv 里全部有cis工具变量的aptamer（deCODE SomaScan，约1,829个）
结局  IgA肾病（Kiryluk 2023，活检确诊）：欧洲 5,556例、亚洲 4,590例、合并 10,146例、去除Immunochip 8,139例
      FinnGen R12：IgA肾病 822例、糖尿病肾病 5,579例、肾小球肾炎 3,905例、肾病综合征 1,106例、肾小球疾病 7,426例
方法  cis-MR（单位点Wald / 多位点固定效应IVW），每个结局内Bonferroni与FDR；
      显著或靠前的蛋白再抽cis区全部位点做共定位（先验与前面各步相同）
主分析用欧洲IgA肾病（工具变量来自冰岛人群，LD结构匹配）；亚洲与合并结果用于看方向是否一致。
MHC区（6号染色体25到34 Mb）的蛋白单独标出：IgA肾病在HLA区信号极强，该区LD复杂，cis-MR不可靠。

不联网。两遍扫描9个结局文件，约20分钟。
输出  extracted\\kd_instrument_assoc.tsv ; results\\kidney_disease_protein_mr.tsv ;
      results\\kidney_disease_coloc.tsv ; logs\\kidney_disease_mr_report.txt
"""
import csv
import glob
import gzip
import importlib.util
import math
import os
import sys
import time
import traceback

import numpy as np

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
OLD_PROJ = os.environ.get("MR_OLD_PROJ", r"D:\肥胖国自然")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "kidney_disease_mr_report.txt")
CACHE = os.path.join(EXT, "kd_instrument_assoc.tsv")
CIS_INSTR = os.path.join(OLD_PROJ, "extracted", "cis_instruments.tsv")
CIS_DIR = os.path.join(OLD_PROJ, "extracted", "cis")
NO = os.path.join(PROJ, "raw", "outcome")
K = os.path.join(NO, "igan_kiryluk2023")
OUTFILES = [
    ("IGAN_EUR", "IgA肾病 欧洲 5,556例", os.path.join(K, "IgAN_Combined_metaanalysis_European_only.txt"), "kiryluk"),
    ("IGAN_ASN", "IgA肾病 亚洲 4,590例", os.path.join(K, "IgAN_Combined_metaanalysis_Asian_only.txt"), "kiryluk"),
    ("IGAN_ALL", "IgA肾病 合并 10,146例", os.path.join(K, "IgAN_Combined_metaanalysis.txt"), "kiryluk"),
    ("IGAN_NOIC", "IgA肾病 去Immunochip 8,139例", os.path.join(K, "IgAN_Combined_metaanalysis_without_immunochip.txt"), "kiryluk"),
    ("IGAN_FG", "IgA肾病 FinnGen 822例", os.path.join(NO, "finngen_r12_n14_iga_nephropathy", "finngen_R12_N14_IGA_NEPHROPATHY.gz"), "finngen"),
    ("DKD_FG", "糖尿病肾病 FinnGen 5,579例", os.path.join(NO, "finngen_r12_dm_nephropathy", "finngen_R12_DM_NEPHROPATHY.gz"), "finngen"),
    ("GN_FG", "肾小球肾炎 FinnGen 3,905例", os.path.join(NO, "finngen_r12_n14_glomer_nephritis", "finngen_R12_N14_GLOMER_NEPHRITIS.gz"), "finngen"),
    ("NS_FG", "肾病综合征 FinnGen 1,106例", os.path.join(NO, "finngen_r12_n14_nephroticsynd", "finngen_R12_N14_NEPHROTICSYND.gz"), "finngen"),
    ("GLOM_FG", "肾小球疾病 FinnGen 7,426例", os.path.join(NO, "finngen_r12_n14_glomerular", "finngen_R12_N14_GLOMERULAR.gz"), "finngen"),
]
LABEL = {o: l for o, l, _p, _f in OUTFILES}
PRIMARY = "IGAN_EUR"
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
MAX_COLOC = 80
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def opener(path):
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if gz else \
        open(path, "r", encoding="utf-8", errors="replace")


def scan(path, fmt, want, first=False):
    d, n = {}, 0
    with opener(path) as f:
        head = f.readline().split()
        low = [h.lower().lstrip("#") for h in head]
        if first:
            log(f"    表头：{head[:14]}")
        try:
            if fmt == "kiryluk":
                irs, ia1, ia2, ib, ise = (low.index(x) for x in ("snp", "a1", "a2", "beta", "se"))
            else:
                irs, ia1, ia2, ib, ise = (low.index(x) for x in ("rsids", "alt", "ref", "beta", "sebeta"))
        except ValueError:
            log(f"    !! 列名对不上：{head}")
            return None
        mx = max(irs, ia1, ia2, ib, ise)
        ic = low.index("chr") if "chr" in low else None
        ip = next((low.index(x) for x in ("bp_hg19", "bp", "pos") if x in low), None)
        use_pos = None            # SNP列不是rsID时，改用 染色体:GRCh37坐标 经1000G面板对到rsID
        p2r = {}
        for line in f:
            n += 1
            v = line.split("\t") if "\t" in line else line.split()
            if len(v) <= mx:
                continue
            if use_pos is None:
                use_pos = fmt == "kiryluk" and not v[irs].startswith("rs") and ic is not None and ip is not None
                if first:
                    log(f"    样例行：{v[:8]}" + ("  → SNP列不是rsID，改按坐标匹配" if use_pos else ""))
                if use_pos:
                    p2r = pos_lookup(want)
            if use_pos:
                rs = p2r.get((v[ic].replace("chr", ""), v[ip]))
                if rs is None:
                    continue
            else:
                rs = v[irs]
            if rs not in want:
                if "," not in rs:
                    continue
                hit = [x for x in rs.split(",") if x in want]
                if not hit:
                    continue
                rs = hit[0]
            try:
                b, se = float(v[ib]), float(v[ise])
            except ValueError:
                continue
            if se > 0 and math.isfinite(b) and rs not in d:
                d[rs] = (v[ia1].strip().upper(), v[ia2].strip().upper(), b, se)
    if first:
        log(f"    扫描 {n:,} 行，抽到 {len(d):,} 个位点")
    return d


_BIM = {}


def pos_lookup(want):
    """{(chr, pos37): rsid}，只含want里的rsID；1000G EUR面板的.bim只读一次"""
    if not _BIM:
        bim = os.path.join(OLD_PROJ, "raw", "reference", "ld", "EUR.bim")
        with open(bim, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                v = line.split()
                if len(v) >= 4:
                    _BIM[v[1]] = (v[0].replace("chr", ""), v[3])
        log(f"    读入1000G面板 {len(_BIM):,} 个位点用于坐标匹配")
    return {_BIM[rs]: rs for rs in want if rs in _BIM}


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


def cis_mr(snps, od):
    num = den = 0.0
    used = []
    for rs, ea, oa, bx, _sx in snps:
        o = od.get(rs)
        if o is None or bx == 0 or {ea, oa} in ({"A", "T"}, {"C", "G"}):
            continue
        s = align(ea, oa, o[0], o[1])
        if s == 0:
            continue
        w = 1.0 / o[3] ** 2
        num += w * bx * s * o[2]
        den += w * bx * bx
        used.append(rs)
    if not used or den <= 0:
        return None
    b, se = num / den, math.sqrt(1.0 / den)
    return b, se, pz(b / se), used


def bh(ps):
    p = np.asarray(ps, float)
    n = len(p)
    o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(adj, 0, 1)
    return out


def run():
    m132 = load("_mod132", "132_gsta1_dm_olink.py")
    cis, meta = {}, {}
    with open(CIS_INSTR, "r", encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames or []
        cc = next((c for c in cols if c.lower() in ("chrom", "chr", "chromosome")), None)
        pc = next((c for c in cols if c.lower() in ("pos", "position")), None)
        for r in rd:
            try:
                b, se = float(r["Beta"]), float(r["SE"])
            except (ValueError, KeyError, TypeError):
                continue
            if se <= 0:
                continue
            key = (r["SeqId"], r["Gene"])
            cis.setdefault(key, []).append((r["rsid"], r["EA"].upper(), r["OA"].upper(), b, se))
            if key not in meta:
                try:
                    meta[key] = (str(r[cc]).replace("chr", "") if cc else "", int(float(r[pc])) if pc else 0)
                except (ValueError, TypeError):
                    meta[key] = ("", 0)
    want = {s[0] for v in cis.values() for s in v}
    log(f"有cis工具变量的aptamer {len(cis):,} 个，工具变量rsID {len(want):,} 个；cis_instruments列：{cols}")

    def is_mhc(key):
        c, p = meta.get(key, ("", 0))
        return c == "6" and 25_000_000 <= p <= 34_500_000

    # ---- A 第一遍：抽工具变量
    log("")
    log("A. 从结局文件抽工具变量")
    assoc = {}
    if os.path.exists(CACHE):
        with open(CACHE, "r", encoding="utf-8") as f:
            f.readline()
            for line in f:
                o, rs, ea, oa, b, se = line.rstrip("\n").split("\t")
                assoc.setdefault(o, {})[rs] = (ea, oa, float(b), float(se))
    for o, lab, path, fmt in OUTFILES:
        if o in assoc:
            log(f"  {o:<10} 已有缓存 {len(assoc[o]):,} 个位点")
            continue
        log(f"  {o}（{lab}）← {path}")
        if not os.path.exists(path):
            log("    !! 文件不存在，跳过")
            continue
        t0 = time.time()
        d = scan(path, fmt, want, first=True)
        if d is None:
            continue
        assoc[o] = d
        log(f"    用时 {(time.time() - t0) / 60:.1f} 分钟")
        with open(CACHE, "w", encoding="utf-8", newline="") as f:
            f.write("outcome\trsid\tEA\tOA\tbeta\tse\n")
            for oo, dd in assoc.items():
                for rs, (ea, oa, b, se) in dd.items():
                    f.write(f"{oo}\t{rs}\t{ea}\t{oa}\t{b:.8g}\t{se:.8g}\n")

    # ---- B cis-MR
    log("")
    log("B. cis-MR：每SD蛋白对应的ln(OR)")
    res = {}
    for o, lab, _p, _f in OUTFILES:
        d = assoc.get(o)
        if not d:
            continue
        rows = []
        for key, snps in cis.items():
            r = cis_mr(snps, d)
            if r:
                rows.append(dict(outcome=o, SeqId=key[0], Gene=key[1], nSNP=len(r[3]), beta=r[0], se=r[1], p=r[2],
                                 rsids=";".join(r[3]), MHC=int(is_mhc(key))))
        if not rows:
            continue
        q = bh([r["p"] for r in rows])
        for r, qv in zip(rows, q):
            r["FDR"] = float(qv)
            r["bonferroni"] = int(r["p"] < 0.05 / len(rows))
            res[(o, r["SeqId"])] = r
        rows.sort(key=lambda r: r["p"])
        nb = [r for r in rows if r["bonferroni"]]
        nf = [r for r in rows if r["FDR"] < 0.05]
        log(f"  {o:<10} {lab}：可检验 {len(rows):,}；Bonferroni {len(nb)}（其中MHC区 {sum(r['MHC'] for r in nb)}）；FDR<0.05 {len(nf)}")
        for r in rows[:15]:
            log(f"      {r['Gene']:<10} {r['SeqId']:<10} OR={math.exp(r['beta']):.2f}  beta={r['beta']:+.3f}  P={r['p']:.1e}  "
                f"FDR={r['FDR']:.1e}  nSNP={r['nSNP']}{'  [MHC]' if r['MHC'] else ''}")
    allrows = sorted(res.values(), key=lambda r: (r["outcome"], r["p"]))
    with open(os.path.join(RES, "kidney_disease_protein_mr.tsv"), "w", encoding="utf-8", newline="") as f:
        cols2 = ["outcome", "SeqId", "Gene", "nSNP", "beta", "se", "p", "FDR", "bonferroni", "MHC", "rsids"]
        f.write("\t".join(cols2) + "\n")
        for r in allrows:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols2) + "\n")

    # ---- C IgA肾病跨人群一致性
    log("")
    log("C. IgA肾病候选：欧洲FDR<0.05或合并FDR<0.05的蛋白，在各数据集中的效应（MHC区已标）")
    igsets = ["IGAN_EUR", "IGAN_ASN", "IGAN_ALL", "IGAN_NOIC", "IGAN_FG"]
    seqs = sorted({s for (o, s), r in res.items() if o in ("IGAN_EUR", "IGAN_ALL") and r["FDR"] < 0.05},
                  key=lambda s: min(res[(o, s)]["p"] for o in ("IGAN_EUR", "IGAN_ALL") if (o, s) in res))
    log("    基因       aptamer     " + "  ".join(f"{o:<16}" for o in igsets) + " 方向")
    for s in seqs[:60]:
        any_r = next(res[(o, s)] for o in igsets if (o, s) in res)
        cells, signs = [], []
        for o in igsets:
            r = res.get((o, s))
            if r:
                cells.append(f"{r['beta']:+.2f}({r['p']:.0e})".ljust(16))
                if o in ("IGAN_EUR", "IGAN_ASN", "IGAN_FG"):
                    signs.append(1 if r["beta"] > 0 else -1)
            else:
                cells.append("-".ljust(16))
        cons = "一致" if len(signs) >= 2 and len(set(signs)) == 1 else ("不一致" if len(signs) >= 2 else "")
        log(f"    {any_r['Gene']:<10} {s:<11} " + "  ".join(cells) + f" {cons}{'  [MHC]' if any_r['MHC'] else ''}")
    if not seqs:
        log("    没有FDR<0.05的蛋白")

    # ---- D 共定位
    log("")
    log("D. 共定位（cis区全部位点 × 各结局）")
    pick = []
    for o, _l, _p, _f in OUTFILES:
        rs_o = sorted([r for (oo, _s), r in res.items() if oo == o], key=lambda r: r["p"])
        pick += [(r["SeqId"], r["Gene"]) for r in rs_o if r["bonferroni"] or r["FDR"] < 0.05]
        if o in ("IGAN_EUR", "IGAN_ALL"):
            pick += [(r["SeqId"], r["Gene"]) for r in rs_o[:15]]
    seen, targets = set(), []
    for k in pick:
        if k not in seen:
            seen.add(k)
            targets.append(k)
    targets = targets[:MAX_COLOC]
    log(f"  进入共定位的aptamer {len(targets)} 个")
    cisreg, missing = {}, []
    for seq, gene in targets:
        p = os.path.join(CIS_DIR, f"{seq}_{gene}.tsv.gz")
        if not os.path.exists(p):
            missing.append(f"{gene}({seq})")
            continue
        for s2, d in m132.read_cis_all(gene):
            if s2 == seq:
                cisreg[(seq, gene)] = d
    if missing:
        log(f"  !! D盘没有cis区文件、无法做共定位的 {len(missing)} 个：{'，'.join(missing[:40])}")
    want2 = set()
    for d in cisreg.values():
        want2 |= set(d)
    crow = []
    if want2:
        log(f"  第二遍扫描：cis区位点 {len(want2):,} 个")
        for o, lab, path, fmt in OUTFILES:
            if not os.path.exists(path):
                continue
            od = scan(path, fmt, want2)
            if od is None:
                continue
            for (seq, gene), d in cisreg.items():
                b1, s1, b2, s2 = m132.harmonised(d, od)
                if len(b1) < 50:
                    continue
                pp = m132.coloc(b1, s1, b2, s2)
                r = res.get((o, seq))
                crow.append(dict(outcome=o, SeqId=seq, Gene=gene, n=len(b1), PP3=float(pp[3]), PP4=float(pp[4]),
                                 mr_beta=r["beta"] if r else float("nan"), mr_p=r["p"] if r else float("nan"),
                                 MHC=int(is_mhc((seq, gene)))))
            log(f"    {o} 完成")
        with open(os.path.join(RES, "kidney_disease_coloc.tsv"), "w", encoding="utf-8", newline="") as f:
            c3 = ["outcome", "SeqId", "Gene", "n", "PP3", "PP4", "mr_beta", "mr_p", "MHC"]
            f.write("\t".join(c3) + "\n")
            for r in crow:
                f.write("\t".join(f"{r[c]:.5g}" if isinstance(r[c], float) else str(r[c]) for c in c3) + "\n")
        log("  PP4>=0.5 的组合（按PP4排序）：")
        good = sorted([r for r in crow if r["PP4"] >= 0.5], key=lambda r: -r["PP4"])
        for r in good[:60]:
            log(f"    {r['Gene']:<10} {r['SeqId']:<10} {r['outcome']:<10} PP4={r['PP4']:.2f} PP3={r['PP3']:.2f}  "
                f"MR beta={r['mr_beta']:+.3f} P={r['mr_p']:.1e}{'  [MHC]' if r['MHC'] else ''}")
        if not good:
            log("    没有")
        log("  IgA肾病：MR显著（欧洲或合并FDR<0.05）蛋白的PP4一览")
        byk = {}
        for r in crow:
            byk.setdefault((r["SeqId"], r["Gene"]), {})[r["outcome"]] = r["PP4"]
        for s in seqs[:60]:
            k = next((kk for kk in byk if kk[0] == s), None)
            if k:
                log(f"    {k[1]:<10} {s:<11} " + "  ".join(f"{o}:{byk[k].get(o, float('nan')):.2f}" for o in igsets))


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    for d in (LOG_DIR, RES, EXT):
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
