# -*- coding: utf-8 -*-
"""149步：EIF2AK3(PERK) 与 ST6GAL1 的工具变量核查
A. 领头pQTL与错义变异/次级信号的LD（1000G欧洲，直接读bed）
B. 多个独立cis工具变量逐个看：IgA肾病（4份）与血清IgA水平是否同向
C. GTEx：领头pQTL是否也是该基因的eQTL（需要联网；失败不影响A、B）
只读 F:\\孟德尔 与 D:\\肥胖国自然\\raw\\reference\\ld
"""
import csv
import importlib.util
import json
import math
import os
import sys
import urllib.parse

import numpy as np

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
OLD = os.environ.get("MR_OLD_PROJ", r"D:\肥胖国自然")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "perk_st6gal1_instruments_report.txt")
INSTR = os.path.join(EXT, "olink_cis_instruments.tsv")
BFILE = os.path.join(OLD, "raw", "reference", "ld", "EUR")
K = os.path.join(PROJ, "raw", "outcome", "igan_kiryluk2023")
OUTCOMES = [("IgAN欧洲", os.path.join(K, "IgAN_Combined_metaanalysis_European_only.txt")),
            ("IgAN亚洲", os.path.join(K, "IgAN_Combined_metaanalysis_Asian_only.txt")),
            ("IgAN合并", os.path.join(K, "IgAN_Combined_metaanalysis.txt")),
            ("IgAN去IC", os.path.join(K, "IgAN_Combined_metaanalysis_without_immunochip.txt")),
            ("IgA水平", os.path.join(PROJ, "raw", "outcome", "serum_iga_liu2022", "META.IGA.LEVELS.ALL.COMBINED.txt"))]
GENES = ["EIF2AK3", "ST6GAL1"]
# PERK单倍型B的三个错义变异（S136C、R166Q、S704A）
EXTRA = {"EIF2AK3": ["rs867529", "rs13045", "rs1805165"], "ST6GAL1": []}
TISSUES = ["Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Spleen", "Small_Intestine_Terminal_Ileum",
           "Liver", "Kidney_Cortex"]
_out, _rows = [], []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2))


def read_instr():
    d = {}
    with open(INSTR, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Gene"] in GENES:
                d.setdefault(r["Gene"], []).append(
                    dict(rsid=r["rsid"], chr=r["CHR"], pos37=int(r["POS37"]), pos38=r["POS38"], oa=r["OA"].upper(),
                         ea=r["EA"].upper(), b=float(r["Beta"]), se=float(r["SE"]), src=r["source"]))
    return d


class Bed:
    """直接读plink bed，取个别位点的基因型"""

    def __init__(self, bfile, want_rs, want_pos):
        self.ok = all(os.path.exists(bfile + e) for e in (".bed", ".bim", ".fam"))
        self.idx = {}
        if not self.ok:
            return
        with open(bfile + ".fam", "r") as f:
            self.n = sum(1 for _ in f)
        self.bps = (self.n + 3) // 4
        with open(bfile + ".bim", "r") as f:
            for i, line in enumerate(f):
                v = line.split()
                if v[1] in want_rs or (v[0], v[3]) in want_pos:
                    key = v[1] if v[1] in want_rs else want_pos[(v[0], v[3])]
                    self.idx.setdefault(key, (i, v[4].upper(), v[5].upper()))
        self.fh = open(bfile + ".bed", "rb")
        if self.fh.read(3) != b"\x6c\x1b\x01":
            self.ok = False

    def geno(self, key):
        """返回 (A1剂量数组含nan, A1, A2)"""
        if not self.ok or key not in self.idx:
            return None
        i, a1, a2 = self.idx[key]
        self.fh.seek(3 + i * self.bps)
        raw = np.frombuffer(self.fh.read(self.bps), dtype=np.uint8)
        two = np.stack([(raw >> s) & 3 for s in (0, 2, 4, 6)], axis=1).reshape(-1)[:self.n]
        g = np.array([2.0, np.nan, 1.0, 0.0])[two]
        return g, a1, a2

    def r(self, k1, k2):
        x, y = self.geno(k1), self.geno(k2)
        if x is None or y is None:
            return None
        ok = ~np.isnan(x[0]) & ~np.isnan(y[0])
        if ok.sum() < 50 or x[0][ok].std() == 0 or y[0][ok].std() == 0:
            return None
        return float(np.corrcoef(x[0][ok], y[0][ok])[0, 1]), x[1], y[1]


def part_a(instr, bed):
    log("A. 与领头pQTL的LD（1000G欧洲）")
    if not bed.ok:
        log(f"  找不到 {BFILE}.bed/.bim/.fam，跳过LD")
        return
    for g in GENES:
        lead = instr[g][0]
        log(f"  {g} 领头 {lead['rsid']}（效应等位基因 {lead['ea']}，蛋白beta {lead['b']:+.3f}）")
        for rs in [i["rsid"] for i in instr[g][1:]] + EXTRA[g]:
            x = bed.r(lead["rsid"], rs)
            tag = "错义" if rs in EXTRA[g] else "次级信号"
            if x is None:
                log(f"    {rs:<12} {tag}  参考面板里没有")
            else:
                log(f"    {rs:<12} {tag}  r={x[0]:+.3f}  r2={x[0] ** 2:.3f}  （{lead['rsid']}-{x[1]} 与 {rs}-{x[2]} 的相关）")


def part_b(instr, bed, m147, m98):
    log("")
    log("B. 逐个工具变量：Wald = 每SD蛋白对应的结局效应（IgAN为log OR，IgA水平为标准化单位）")
    regions = {}
    for g in GENES:
        for i in instr[g]:
            regions.setdefault(i["chr"], []).append((i["pos37"], i["pos37"]))
    for name, path in OUTCOMES:
        if not os.path.exists(path):
            log(f"  {name}：找不到文件")
            continue
        d, n = m147.scan_igan(path, regions)
        log(f"  {name}（扫描 {n:,} 行）")
        for g in GENES:
            keep = []
            for i in instr[g]:
                o = d.get((i["chr"], i["pos37"]))
                if o is None:
                    log(f"    {g:<8} {i['rsid']:<12} 结局里没有")
                    continue
                s = m147.sgn(i["ea"], i["oa"], o[0], o[1])
                if s == 0:
                    log(f"    {g:<8} {i['rsid']:<12} 等位基因对不上")
                    continue
                by, bse = s * o[2], o[3]
                w, wse = by / i["b"], bse / abs(i["b"])
                log(f"    {g:<8} {i['rsid']:<12} {i['src']:<10} 蛋白beta={i['b']:+.3f}  Wald={w:+.3f} (se {wse:.3f}) P={pz(w / wse):.2e}")
                _rows.append([g, name, i["rsid"], i["src"], f"{i['b']:.5g}", f"{w:.5g}", f"{wse:.5g}", f"{pz(w / wse):.4g}"])
                keep.append((abs(i["b"] / i["se"]), i, by, bse))
            # 按蛋白z从大到小，r2<0.1才保留
            keep.sort(key=lambda t: -t[0])
            sel = []
            for t in keep:
                ok = True
                for u in sel:
                    x = bed.r(t[1]["rsid"], u[1]["rsid"]) if bed.ok else None
                    if x is None or x[0] ** 2 >= 0.1:
                        ok = False
                        break
                if ok:
                    sel.append(t)
            if len(sel) >= 2:
                bx = np.array([t[1]["b"] for t in sel]); by = np.array([t[2] for t in sel]); bse = np.array([t[3] for t in sel])
                b, se, p, _q, qp, i2 = m98.ivw(bx, by, bse)
                log(f"    {g:<8} IVW（r2<0.1的 {len(sel)} 个：{','.join(t[1]['rsid'] for t in sel)}）={b:+.3f} (se {se:.3f}) P={p:.2e}  异质性P={qp:.2g}")
                _rows.append([g, name, "IVW_%d" % len(sel), "pruned_r2_0.1", "", f"{b:.5g}", f"{se:.5g}", f"{p:.4g}"])
            else:
                log(f"    {g:<8} 独立工具变量不足2个（或缺LD），不做IVW")


def part_c(instr, m141):
    log("")
    log("C. GTEx v8：领头pQTL × 本基因表达")
    for g in GENES:
        lead = instr[g][0]
        raw = m141.http("https://gtexportal.org/api/v2/reference/gene?" + urllib.parse.urlencode(
            {"geneId": g, "gencodeVersion": "v26", "genomeBuild": "GRCh38/hg38"}))
        if raw is None:
            log(f"  {g}：查不到gencodeId（网络？）")
            continue
        try:
            gid = json.loads(raw.decode("utf-8", "replace"))["data"][0]["gencodeId"]
        except (ValueError, KeyError, IndexError):
            log(f"  {g}：gencodeId解析失败")
            continue
        for tis in TISSUES:
            got = False
            for ref, alt in ((lead["oa"], lead["ea"]), (lead["ea"], lead["oa"])):
                vid = f"chr{lead['chr']}_{lead['pos38']}_{ref}_{alt}_b38"
                raw = m141.http("https://gtexportal.org/api/v2/association/dyneqtl?" + urllib.parse.urlencode(
                    {"gencodeId": gid, "variantId": vid, "tissueSiteDetailId": tis, "datasetId": "gtex_v8"}), tries=1)
                if raw is None:
                    continue
                try:
                    o = json.loads(raw.decode("utf-8", "replace"))
                except ValueError:
                    continue
                if o.get("nes") is None:
                    continue
                log(f"  {g:<8} {lead['rsid']} 在 {tis:<36} NES={float(o['nes']):+.3f} P={float(o['pValue']):.2e}  （NES对应 {alt}；蛋白效应等位基因 {lead['ea']}，蛋白beta {lead['b']:+.3f}）")
                got = True
                break
            if not got:
                log(f"  {g:<8} {lead['rsid']} 在 {tis:<36} 没取到")


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    m147 = load("_m147", "147_olink_igan_coloc.py")
    m98 = load("_m98", "98_mr_step1_panukb.py")
    instr = read_instr()
    for g in GENES:
        if g not in instr:
            log(f"工具变量表里没有 {g}")
            return
        instr[g].sort(key=lambda i: i["src"] != "ST9_lead")
    want = set(r for g in GENES for r in [i["rsid"] for i in instr[g]] + EXTRA[g])
    bed = Bed(BFILE, want, {})
    for fn, tag in ((lambda: part_a(instr, bed), "A"), (lambda: part_b(instr, bed, m147, m98), "B")):
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            import traceback
            log(f"!! {tag}部分异常：{type(e).__name__}: {e}")
            log(traceback.format_exc())
    if os.environ.get("MR_SKIP_NET"):
        log("\nC. 跳过联网部分")
    else:
        try:
            m141 = load("_m141", "141_gsta1_annotation_probe.py")
            m141.log = log
            part_c(instr, m141)
        except Exception as e:  # noqa: BLE001
            log(f"!! C部分异常：{type(e).__name__}: {e}")
    with open(os.path.join(RES, "perk_st6gal1_instruments.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Gene", "outcome", "instrument", "source", "protein_beta", "wald", "se", "p"])
        w.writerows(_rows)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
