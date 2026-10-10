# -*- coding: utf-8 -*-
"""152步：ST6GAL1位点的两个信号拆开看
  S1 = 循环ST6GAL1的领头pQTL  chr3:186,739,677（rs4686837）
  S2 = IgG唾液酸化的最强位点    chr3:186,727,170
对每个性状给出两个位点各自的边际z，以及互相校正后的z（用1000G欧洲LD，只对欧洲为主的数据做校正）。
只读，不联网，约5分钟（用151步的缓存）。
"""
import csv
import glob
import importlib.util
import math
import os
import sys
import traceback

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PROJ, "logs")
RES = os.path.join(PROJ, "results")
REPORT = os.path.join(LOG_DIR, "st6gal1_two_signals_report.txt")
CACHE = os.path.join(PROJ, "extracted", "igg_glycome_st6gal1_window.tsv")
IGA = os.path.join(PROJ, "raw", "outcome", "serum_iga_liu2022", "META.IGA.LEVELS.ALL.COMBINED.txt")
CHROM, S1, S2 = "3", 186739677, 186727170
TRAITS = [("IGP29", "FA2G2S1占比（岩藻糖化双半乳糖结构的单唾液酸化）"), ("IGP24", "岩藻糖化半乳糖化结构的唾液酸化率"),
          ("IGP15", "GP16 FA2G1S1"), ("IGP13", "GP14 FA2G2（未唾液酸化底物）"), ("IGP14", "GP15 FA2BG2"),
          ("IGP54", "GP15/中性糖链（FA2BG2）"), ("IGP52", "GP13/中性糖链（A2BG2）"), ("IGP43", "GP4/中性糖链（FA2，无半乳糖）"),
          ("IGP55", "中性糖链中无半乳糖结构合计"), ("IGP57", "中性糖链中双半乳糖结构合计")]
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


def run():
    m = load("_m147", "147_olink_igan_coloc.py")
    m149 = load("_m149", "149_perk_st6gal1_instruments.py")
    bed = m149.Bed(m149.BFILE, set(), {(CHROM, str(S1)): "S1", (CHROM, str(S2)): "S2"})
    r, ref = None, {}
    if bed.ok and "S1" in bed.idx and "S2" in bed.idx:
        x = bed.r("S1", "S2")
        ref = {"S1": bed.idx["S1"][1:], "S2": bed.idx["S2"][1:]}
        if x:
            r = x[0]
            log(f"LD（1000G欧洲）：S1-{x[1]} 与 S2-{x[2]} 的 r={r:+.3f}，r2={r * r:.3f}")
            for k in ("S1", "S2"):
                g = bed.geno(k)[0]
                import numpy as np
                log(f"  {k} 参考等位基因 {ref[k][0]}（另一 {ref[k][1]}），{ref[k][0]} 的频率 {np.nanmean(g) / 2:.3f}")
    if r is None:
        log("!! 参考面板里找不到这两个位点，无法做互相校正，只给边际结果")

    def line(name, d, cond, note=""):
        """d: {pos: (效应等位基因, 另一等位基因, beta, se)}"""
        z = {}
        for k, p in (("S1", S1), ("S2", S2)):
            o = d.get(p)
            if o is None:
                continue
            a = ref.get(k)
            s = m.sgn(o[0], o[1], a[0], a[1]) if a else 1
            if s == 0:
                continue
            z[k] = (s * o[2] / o[3], s * o[2], o[3])
        t = f"  {name:<16}"
        for k in ("S1", "S2"):
            t += f"  {k}: beta={z[k][1]:+.4f} z={z[k][0]:+6.1f} P={pz(z[k][0]):.1e}" if k in z else f"  {k}: 没有"
        c1 = c2 = float("nan")
        if cond and r is not None and len(z) == 2:
            d_ = math.sqrt(1 - r * r)
            c1, c2 = (z["S1"][0] - r * z["S2"][0]) / d_, (z["S2"][0] - r * z["S1"][0]) / d_
            t += f"   | 互相校正后 S1 z={c1:+6.1f}（P={pz(c1):.1e}）  S2 z={c2:+6.1f}（P={pz(c2):.1e}）"
        log(t + ("   " + note if note else ""))
        _rows.append([name] + [f"{z[k][i]:.5g}" if k in z else "" for k in ("S1", "S2") for i in (1, 2, 0)] + [f"{c1:.4g}", f"{c2:.4g}"])

    log("")
    log("（beta、z的方向都对应上面列出的参考等位基因）")
    log("A. 循环ST6GAL1（Olink，欧洲）")
    tars = [t for t in glob.glob(os.path.join(m.TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith("ST6GAL1_")]
    if tars:
        ol = m.read_olink(tars[0], CHROM, min(S1, S2) - 10, max(S1, S2) + 10)
        line("ST6GAL1蛋白", {p: (v[1], v[0], v[2], v[3]) for p, v in ol.items()}, True)
    log("")
    log("B. IgG糖链性状（欧洲）")
    gly = {}
    if os.path.exists(CACHE):
        with open(CACHE, "r", encoding="utf-8") as f:
            for q in csv.DictReader(f, delimiter="\t"):
                if int(q["pos37"]) in (S1, S2):
                    gly.setdefault(q["trait"], {})[int(q["pos37"])] = (q["ea"], q["oa"], float(q["beta"]), float(q["se"]))
    else:
        log("  没有151步的缓存，先跑151步")
    for t, desc in TRAITS:
        if t in gly:
            line(t, gly[t], True, desc)
    log("")
    log("C. IgA肾病与血清IgA水平（亚洲、合并数据的LD不同，不做校正）")
    for name, path in list(m.IGAN) + [("IgA水平", IGA)]:
        if not os.path.exists(path):
            continue
        d, _n = m.scan_igan(path, {CHROM: [(S1, S1), (S2, S2)]})
        line("IgAN " + name if name != "IgA水平" else name, {k[1]: v for k, v in d.items()}, name in ("欧洲", "IgA水平"))
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "st6gal1_two_signals.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", "S1_beta", "S1_se", "S1_z", "S2_beta", "S2_se", "S2_z", "S1_z_cond_on_S2", "S2_z_cond_on_S1"])
        w.writerows(_rows)


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    try:
        run()
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
