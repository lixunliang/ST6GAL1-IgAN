# -*- coding: utf-8 -*-
"""156步：用东亚参考人群的LD重做ST6GAL1两个信号的互相校正
  S1 rs4686837（chr3:186,739,677） S2 rs7621161（chr3:186,727,170）
需要 D:\\肥胖国自然\\raw\\reference\\ld\\EAS.bed/.bim/.fam（与EUR同一套参考数据）。没有的话报告里会说明。
只读，不联网，约3分钟。
"""
import glob
import importlib.util
import math
import os
import sys
import traceback

import numpy as np

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
OLD = os.environ.get("MR_OLD_PROJ", r"D:\肥胖国自然")
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(PROJ, "logs", "st6gal1_eas_ld_report.txt")
LD_DIR = os.path.join(OLD, "raw", "reference", "ld")
CHROM, S1, S2 = "3", 186739677, 186727170
# 155步：BBJ（日本）里两位点A等位基因的效应（beta, P）
BBJ = [("非白蛋白蛋白（球蛋白）", (-0.026, 1.6e-10), (-0.014, 9.6e-4)),
       ("白蛋白/球蛋白比", (0.023, 6.8e-9), (0.013, 2.6e-3)),
       ("2型糖尿病（多人种）", (-0.013, 3.6e-5), (-0.0092, 2.8e-3)),
       ("单核细胞计数", (-0.015, 9.1e-4), (-0.012, 1.2e-2)),
       ("白细胞计数", (-0.0083, 1.8e-2), (-0.013, 4.7e-4))]
_out = []


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


def zfromp(beta, p):
    from statistics import NormalDist
    return math.copysign(abs(NormalDist().inv_cdf(p / 2)), beta)


def cond(z1, z2, r):
    d = math.sqrt(1 - r * r)
    return (z1 - r * z2) / d, (z2 - r * z1) / d


def run():
    m = load("_m147", "147_olink_igan_coloc.py")
    m149 = load("_m149", "149_perk_st6gal1_instruments.py")
    log(f"参考数据文件夹 {LD_DIR} 里的bed文件：" + ", ".join(sorted(os.path.basename(p) for p in glob.glob(os.path.join(LD_DIR, "*.bed")))))
    rs, refs = {}, {}
    for pop in ("EUR", "EAS"):
        bed = m149.Bed(os.path.join(LD_DIR, pop), set(), {(CHROM, str(S1)): "S1", (CHROM, str(S2)): "S2"})
        if not bed.ok:
            log(f"  {pop}：没有 {pop}.bed/.bim/.fam")
            continue
        x = bed.r("S1", "S2")
        if not x:
            log(f"  {pop}：参考面板里缺这两个位点之一（{sorted(bed.idx)}）")
            continue
        rs[pop] = x[0]
        refs[pop] = {k: bed.idx[k][1:] for k in ("S1", "S2")}
        f1, f2 = (np.nanmean(bed.geno(k)[0]) / 2 for k in ("S1", "S2"))
        log(f"  {pop}（{bed.n}人）：S1-{x[1]} 与 S2-{x[2]} 的 r={x[0]:+.3f}，r2={x[0] ** 2:.3f}；频率 S1-{x[1]} {f1:.3f}，S2-{x[2]} {f2:.3f}")
    if "EAS" not in rs:
        log("")
        log("没有东亚LD，无法做东亚校正。需要IEU的1000G参考包里的 EAS.bed/.bim/.fam 放到上面的文件夹。")
        return
    log("")
    log("A. IgA肾病（方向对应各人群参考面板的第一个等位基因，见上）")
    for name, path in m.IGAN:
        if not os.path.exists(path):
            continue
        d, _n = m.scan_igan(path, {CHROM: [(S1, S1), (S2, S2)]})
        for pop in (["EAS"] if name == "亚洲" else ["EUR"] if name == "欧洲" else ["EUR", "EAS"]):
            if pop not in rs:
                continue
            z = {}
            for k, p in (("S1", S1), ("S2", S2)):
                o = d.get((CHROM, p))
                if o is None:
                    continue
                s = m.sgn(o[0], o[1], refs[pop][k][0], refs[pop][k][1])
                if s:
                    z[k] = s * o[2] / o[3]
            if len(z) < 2:
                log(f"  IgAN {name}：缺位点")
                continue
            c1, c2 = cond(z["S1"], z["S2"], rs[pop])
            log(f"  IgAN {name:<4} 用{pop} LD：边际 S1 z={z['S1']:+.1f} S2 z={z['S2']:+.1f}；校正后 S1 z={c1:+.1f}（P={pz(c1):.1e}） S2 z={c2:+.1f}（P={pz(c2):.1e}）")
    log("")
    log("B. 日本生物样本库（155步的数值，A等位基因；用EAS LD校正）")
    r = rs["EAS"] * (1 if (refs["EAS"]["S1"][0] == "A") == (refs["EAS"]["S2"][0] == "A") else -1)
    for name, (b1, p1), (b2, p2) in BBJ:
        z1, z2 = zfromp(b1, p1), zfromp(b2, p2)
        c1, c2 = cond(z1, z2, r)
        log(f"  {name:<18} 边际 S1 z={z1:+.1f} S2 z={z2:+.1f}；校正后 S1 z={c1:+.1f}（P={pz(c1):.1e}） S2 z={c2:+.1f}（P={pz(c2):.1e}）")


def main():
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
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
