#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
148_st6gal1_chain.py  —  把"暴露 → 蛋白 → IgA肾病"的链条在ST6GAL1上补全（2026-10-07）

主线蛋白 ST6GAL1；对照 EIF2AK3（PERK，备选）、TNFSF13（APRIL）、TNFRSF13B（TACI）、CFH、ITGAM、IL7R。
A  下游：蛋白 → 血清IgA水平（Liu 2022，41,263人；另有不含deCODE的版本）
   领头cis-pQTL的Wald比值 + ±500 kb共定位。APRIL应明显升高IgA水平（阳性对照）。
B  反方向：IgA肾病易感性 → 蛋白。从IgA肾病合并数据取全基因组显著、相距1 Mb以上的位点作工具变量
   （剔除MHC区和该蛋白基因±1 Mb内的位点），在Olink全基因组数据里做IVW。
   若显著，说明"患病使蛋白改变"，与"蛋白致病"是两回事。
C  上游：5个代谢暴露（空腹血糖、HbA1c、空腹胰岛素、2型糖尿病、收缩压）→ 这批Olink蛋白
   （复用122步的全基因组扫描；146步新下的25个蛋白都会算，约40分钟），并列出SomaScan里的对应结果。
不联网。A、B约15分钟，C约40分钟。
输出  results\\st6gal1_chain.tsv ; logs\\st6gal1_chain_report.txt ; results\\expo_olink_step1.tsv（追加）
"""
import csv
import glob
import importlib.util
import math
import os
import sys
import traceback

import numpy as np

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "st6gal1_chain_report.txt")
IGA_DIR = os.path.join(PROJ, "raw", "outcome", "serum_iga_liu2022")
IGA = [("IgA水平(全部)", os.path.join(IGA_DIR, "META.IGA.LEVELS.ALL.COMBINED.txt")),
       ("IgA水平(不含deCODE)", os.path.join(IGA_DIR, "META.IGA.LEVELS.non-DECODE.txt"))]
IGAN_ALL = os.path.join(PROJ, "raw", "outcome", "igan_kiryluk2023", "IgAN_Combined_metaanalysis.txt")
FOCUS = ["ST6GAL1", "EIF2AK3", "TNFSF13", "TNFRSF13B", "CFH", "ITGAM", "IL7R"]
if os.environ.get("MR_TEST_GENES"):
    FOCUS = os.environ["MR_TEST_GENES"].split(",")
EXPOSURES = ["FG", "HBA1C", "FI", "T2D", "SBP"]
_out, _rows = [], []


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


def part_a(m):
    log("A. 蛋白 → 血清IgA水平（Wald为每SD蛋白对应的IgA水平变化，标准化单位）")
    lead, olink, regions = {}, {}, {}
    with open(m.INSTR, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["source"] == "ST9_lead" and r["Gene"] in FOCUS and r["Gene"] not in lead:
                lead[r["Gene"]] = (r["CHR"], int(r["POS37"]), r["rsid"])
    for g in FOCUS:
        tars = [t for t in glob.glob(os.path.join(m.TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith(g + "_")]
        if g not in lead or not tars:
            log(f"  {g}：缺领头位点或tar，跳过")
            continue
        c, p, _rs = lead[g]
        d = m.read_olink(tars[0], c, p - m.WIN, p + m.WIN)
        if d:
            olink[g] = d
            regions.setdefault(c, []).append((p - m.WIN, p + m.WIN))
    for name, path in IGA:
        if not os.path.exists(path):
            log(f"  {name}：找不到 {path}")
            continue
        try:
            idat, n = m.scan_igan(path, regions)
        except Exception as e:  # noqa: BLE001
            log(f"  {name}：读取失败 {type(e).__name__}: {e}")
            continue
        log(f"  {name}：扫描 {n:,} 行，窗口内 {len(idat):,} 个位点")
        for g in FOCUS:
            if g not in olink:
                continue
            c, p, rs = lead[g]
            pos, b1, s1, b2, s2 = [], [], [], [], []
            for q, (a0, a1, b, se) in olink[g].items():
                o = idat.get((c, q))
                if o is None:
                    continue
                s = m.sgn(a1, a0, o[0], o[1])
                if s == 0:
                    continue
                pos.append(q); b1.append(b); s1.append(se); b2.append(s * o[2]); s2.append(o[3])
            if len(pos) < 50:
                log(f"    {g:<10} 共有位点 {len(pos)} 个，跳过")
                continue
            pos = np.array(pos); b1 = np.array(b1); s1 = np.array(s1); b2 = np.array(b2); s2 = np.array(s2)
            pp = m.coloc(b1, s1, b2, s2)
            z1, z2 = b1 / s1, b2 / s2
            i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
            w, wse = b2[i1] / b1[i1], abs(s2[i1] / b1[i1])
            log(f"    {g:<10} Wald={w:+.3f} (se {wse:.3f}) P={pz(w / wse):.1e}  PP3={pp[3]:.2f} PP4={pp[4]:.2f}  "
                f"窗口内IgA水平最小P={pz(z2[i2]):.1e}（该位点蛋白z={z1[i2]:+.1f}，距蛋白最强位点 {abs(int(pos[i1]) - int(pos[i2])) / 1000:.0f} kb）")
            _rows.append([g, "protein_to_" + name, len(pos), f"{w:.5g}", f"{wse:.5g}", f"{pz(w / wse):.4g}", f"{pp[3]:.4g}", f"{pp[4]:.4g}"])
    return lead


def part_b(m, m98, m122, lead):
    log("")
    log("B. 反方向：IgA肾病易感性 → 蛋白（IVW；工具变量来自IgA肾病合并数据）")
    if not os.path.exists(IGAN_ALL):
        log("  找不到IgA肾病合并文件")
        return
    sig = []
    with open(IGAN_ALL, "r", encoding="utf-8", errors="replace") as f:
        low = [h.lower() for h in f.readline().split()]
        ic, ip, ia1, ia2, ib, ise, ipv = (low.index(x) for x in ("chr", "bp_hg19", "a1", "a2", "beta", "se", "p"))
        for line in f:
            v = line.split()
            try:
                p = float(v[ipv])
                if not p < 5e-8:
                    continue
                c, pos, b, se = v[ic], int(v[ip]), float(v[ib]), float(v[ise])
            except (ValueError, IndexError):
                continue
            a1, a2 = v[ia1].upper(), v[ia2].upper()
            if not c.isdigit() or len(a1) != 1 or len(a2) != 1 or {a1, a2} in ({"A", "T"}, {"C", "G"}):
                continue
            if c == "6" and 25_000_000 <= pos <= 34_500_000:
                continue
            sig.append((p, c, pos, a1, a2, b, se))
    sig.sort()
    inst = []
    for s in sig:
        if all(not (s[1] == t[1] and abs(s[2] - t[2]) < 1_000_000) for t in inst):
            inst.append(s)
    log(f"  IgA肾病全基因组显著位点 {len(sig):,} 个，按1 Mb去重并剔除MHC后的工具变量 {len(inst)} 个")
    if len(inst) < 5:
        log("  工具变量太少，跳过")
        return
    key37 = {(c, str(pos)): f"{c}:{pos}" for _p, c, pos, _a1, _a2, _b, _se in inst}
    info = {f"{c}:{pos}": (c, pos, a1, a2, b, se) for _p, c, pos, a1, a2, b, se in inst}
    for g in FOCUS:
        tars = [t for t in glob.glob(os.path.join(m.TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith(g + "_")]
        if not tars:
            continue
        try:
            rows, _cnt = m122.scan_tar(tars[0], key37, {})
        except Exception as e:  # noqa: BLE001
            log(f"  {g}：读tar失败 {type(e).__name__}: {e}")
            continue
        gc, gp = (lead[g][0], lead[g][1]) if g in lead else ("", 0)
        bx, by, bse, used, dropped = [], [], [], set(), 0
        for rid, ea, oa, b, se in rows:
            if rid in used:
                continue
            c, pos, a1, a2, bi, _sei = info[rid]
            if c == gc and abs(pos - gp) < 1_000_000:
                dropped += 1
                used.add(rid)
                continue
            s = m.sgn(ea, oa, a1, a2)      # 蛋白效应对齐到IgAN的a1
            if s == 0:
                continue
            used.add(rid)
            bx.append(bi); by.append(s * b); bse.append(se)
        if len(bx) < 5:
            log(f"  {g:<10} 可用工具变量 {len(bx)} 个，跳过")
            continue
        bx, by, bse = np.array(bx), np.array(by), np.array(bse)
        b, se, p, _q, qp, i2 = m98.ivw(bx, by, bse)
        log(f"  {g:<10} nSNP={len(bx)}（剔除本基因附近 {dropped} 个）  IVW={b:+.4f} (se {se:.4f}) P={p:.2e}  I2={i2:.0f}%")
        _rows.append([g, "IgAN_liability_to_protein", len(bx), f"{b:.5g}", f"{se:.5g}", f"{p:.4g}", "", ""])


def part_c(m122):
    log("")
    log("C. 上游：5个暴露 → Olink蛋白（详见 logs\\expo_olink_step1_report.txt）")
    m122.SETS = [(k, os.path.join(EXT, f"expo_{k}_instruments.tsv")) for k in EXPOSURES]
    m122.OUT = os.path.join(RES, "expo_olink_step1.tsv")
    m122.REPORT = os.path.join(LOG_DIR, "expo_olink_step1_report.txt")
    m122.FIRST = [g.lower() for g in FOCUS]
    m122._out.clear()
    m122.main()
    soma = {}
    for k in EXPOSURES:
        p = os.path.join(RES, f"expo_step1_{k}.tsv")
        if os.path.exists(p):
            with open(p, "r", encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    soma.setdefault((r["Gene"].upper(), k), []).append(r)
    if not os.path.exists(m122.OUT):
        return
    with open(m122.OUT, "r", encoding="utf-8") as f:
        orows = list(csv.DictReader(f, delimiter="\t"))
    log("  基因        暴露    Olink IVW (P) [加权中位数P]          SomaScan IVW (P) [aptamer]")
    for g in FOCUS:
        for r in orows:
            if r["Gene"].split("_")[0] != g:
                continue
            k = r["instrument_set"]
            ss = soma.get((g, k), [])
            try:
                log(f"  {g:<10} {k:<6} {float(r['beta_IVW']):+.4f} ({float(r['p_IVW']):.1e}) [{float(r['p_WM']):.1e}]   "
                    + "；".join(f"{float(s['beta']):+.4f} ({float(s['p']):.1e}) [{s['SeqId']}]" for s in ss))
                _rows.append([g, "exposure_" + k + "_to_protein_Olink", r["nSNP"], r["beta_IVW"], r["se_IVW"], r["p_IVW"], "", ""])
            except ValueError:
                continue
    log("  其余Olink蛋白里 P<1e-3 的暴露→蛋白：")
    for r in orows:
        try:
            if float(r["p_IVW"]) < 1e-3 and r["Gene"].split("_")[0] not in FOCUS:
                log(f"    {r['Gene']:<10} {r['instrument_set']:<6} {float(r['beta_IVW']):+.4f} ({float(r['p_IVW']):.1e})")
        except ValueError:
            continue


def save():
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "st6gal1_chain.tsv"), "w", encoding="utf-8", newline="") as f:
        f.write("Gene\tanalysis\tnSNP\tbeta\tse\tp\tPP3\tPP4\n")
        for r in _rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(LOG_DIR, exist_ok=True)
    m = load("_m147", "147_olink_igan_coloc.py")
    m98 = load("_m98", "98_mr_step1_panukb.py")
    m122 = load("_m122", "122_olink_step1_egfrcys.py")
    lead = {}
    for name, fn in (("A", lambda: part_a(m)), ("B", lambda: part_b(m, m98, m122, lead)), ("C", lambda: part_c(m122))):
        try:
            r = fn()
            if name == "A" and r:
                lead.update(r)
        except Exception as e:  # noqa: BLE001
            log(f"!! {name}部分异常：{type(e).__name__}: {e}")
            for ln in traceback.format_exc().splitlines():
                log("    " + ln)
        save()


if __name__ == "__main__":
    main()
