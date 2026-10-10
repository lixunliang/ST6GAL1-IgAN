# -*- coding: utf-8 -*-
"""160步：ST6GAL1位点的敏感性分析（回应审稿式意见）
A. 2015年汉族GWAS报告的 rs7634389（Li 2015 Nat Commun，P=7.3e-10，风险等位基因C）
   与 S1（rs4686837）、S2（rs7621161）的连锁（1000G欧洲/东亚），及其在IgAN、蛋白、IgG唾液酸化上的效应
B. 共定位先验敏感性：p12 = 1e-6、5e-6、1e-5、5e-5；先验效应SD 0.1/0.15/0.2
C. 条件共定位：先把整段区域对另一个信号做校正，再共定位
   蛋白|S2 × IgAN|S2（应仍为同一信号）；IGP29|S1 × IgAN|S1（应仍为不同信号）
只读，不联网，约10分钟（用151步的糖组缓存、147步的读取函数、149步的LD读取）。
"""
import csv
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
LD_DIR = os.path.join(OLD, "raw", "reference", "ld")
CACHE = os.path.join(PROJ, "extracted", "igg_glycome_st6gal1_window.tsv")
REPORT = os.path.join(PROJ, "logs", "st6gal1_sensitivity_report.txt")
RES = os.path.join(PROJ, "results", "st6gal1_sensitivity.tsv")
CHROM, S1, S2, LEAD = "3", 186739677, 186727170, 186739677
WIN = 500_000
RS2015 = "rs7634389"
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


def coloc(b1, s1, b2, s2, p1=1e-4, p2=1e-4, p12=1e-5, w=0.15):
    W = w * w

    def labf(b, s):
        r = W / (W + s ** 2)
        return 0.5 * (np.log(1 - r) + r * (b / s) ** 2)

    def lse(x):
        mx = np.max(x)
        return mx + np.log(np.sum(np.exp(x - mx)))
    l1, l2 = labf(b1, s1), labf(b2, s2)
    L1, L2, L4 = lse(l1), lse(l2), lse(l1 + l2)
    L3 = L1 + L2 + math.log1p(-math.exp(min(L4 - L1 - L2, -1e-12)))
    lh = np.array([0.0, math.log(p1) + L1, math.log(p2) + L2, math.log(p1) + math.log(p2) + L3, math.log(p12) + L4])
    return np.exp(lh - lse(lh))


def run():
    m = load("_m147", "147_olink_igan_coloc.py")
    m149 = load("_m149", "149_perk_st6gal1_instruments.py")
    lo, hi = LEAD - WIN, LEAD + WIN
    # ---------- 数据
    tars = [t for t in glob.glob(os.path.join(m.TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith("ST6GAL1_")]
    ol = m.read_olink(tars[0], CHROM, lo, hi) if tars else {}
    prot = {p: (v[1], v[0], v[2], v[3]) for p, v in ol.items()}       # (效应, 另一, b, se)
    ig = {}
    for name, path in m.IGAN:
        if os.path.exists(path):
            d, _n = m.scan_igan(path, {CHROM: [(lo, hi)]})
            ig[name] = {k[1]: v for k, v in d.items()}
    gly = {}
    if os.path.exists(CACHE):
        with open(CACHE, "r", encoding="utf-8") as f:
            for q in csv.DictReader(f, delimiter="\t"):
                if q["trait"] in ("IGP29", "IGP54", "IGP24"):
                    gly.setdefault(q["trait"], {})[int(q["pos37"])] = (q["ea"], q["oa"], float(q["beta"]), float(q["se"]))
    log(f"窗口 chr{CHROM}:{lo:,}-{hi:,}：Olink {len(prot):,}；IgAN " + "，".join(f"{k} {len(v):,}" for k, v in ig.items())
        + "；糖组 " + "，".join(f"{k} {len(v):,}" for k, v in gly.items()))

    # ---------- LD面板：窗口内所有位点
    allpos = set(prot) | set().union(*[set(v) for v in ig.values()]) | set().union(*[set(v) for v in gly.values()])
    want_pos = {(CHROM, str(p)): f"p{p}" for p in allpos}
    beds = {}
    for pop in ("EUR", "EAS"):
        b = m149.Bed(os.path.join(LD_DIR, pop), {RS2015}, want_pos)
        if b.ok:
            beds[pop] = b
            log(f"{pop} 参考面板：窗口内找到 {len(b.idx):,} 个位点（{b.n} 人）")
        else:
            log(f"{pop} 参考面板缺失")
    if not beds:
        return

    def key_of(pos):
        return f"p{pos}"

    # ---------- A. rs7634389
    log("")
    log(f"A. {RS2015}（Li 2015：chr3:186,738,421，风险等位基因C，P=7.3e-10，OR 1.13）")
    p15 = None
    for pop, b in beds.items():
        if RS2015 not in b.idx:
            log(f"  {pop}：参考面板里没有 {RS2015}")
            continue
        i, a1, a2 = b.idx[RS2015]
        g = b.geno(RS2015)[0]
        log(f"  {pop}：{RS2015} 参考等位基因 {a1}（另一 {a2}），{a1} 频率 {np.nanmean(g) / 2:.3f}")
        for tag, pos in (("S1 rs4686837", S1), ("S2 rs7621161", S2)):
            x = b.r(RS2015, key_of(pos))
            if x:
                log(f"    与 {tag}：r={x[0]:+.3f}，r2={x[0] ** 2:.3f}（{RS2015}-{x[1]} 与 {x[2]}）")
                _rows.append(["LD", pop, RS2015, tag, f"{x[0]:.4f}", f"{x[0] ** 2:.4f}", "", "", ""])
        p15 = p15 or None
    # 在bim里找它的坐标
    pos15 = None
    for b in beds.values():
        if RS2015 in b.idx:
            with open(os.path.join(LD_DIR, "EUR.bim") if os.path.exists(os.path.join(LD_DIR, "EUR.bim")) else os.path.join(LD_DIR, "EAS.bim"), "r") as f:
                for line in f:
                    v = line.split()
                    if v[1] == RS2015:
                        pos15 = int(v[3])
                        break
            break
    if pos15:
        log(f"  {RS2015} 在参考面板中的GRCh37坐标 chr{CHROM}:{pos15:,}")
        for name, d in list(ig.items()) + [("Olink ST6GAL1", prot)] + [(k, v) for k, v in gly.items()]:
            o = d.get(pos15)
            if o:
                z = o[2] / o[3]
                log(f"    {name:<14} 效应等位基因 {o[0]}：beta={o[2]:+.4f} se={o[3]:.4f} z={z:+.2f} P={pz(z):.1e}")
                _rows.append(["rs7634389_assoc", name, o[0], f"{o[2]:.5g}", f"{o[3]:.5g}", f"{pz(z):.4g}", "", "", ""])
            else:
                log(f"    {name:<14} 没有该位点")

    # ---------- B. 先验敏感性
    def pairs(d1, d2):
        out = []
        for q, a in d1.items():
            o = d2.get(q)
            if o is None:
                continue
            s = m.sgn(a[0], a[1], o[0], o[1])
            if s:
                out.append((q, a[2], a[3], s * o[2], o[3]))
        if not out:
            return None
        q, b1, s1, b2, s2 = (np.array(x) for x in zip(*out))
        return q, b1, s1, b2, s2
    comps = [("蛋白×IgAN合并", prot, ig.get("合并")), ("蛋白×IgAN亚洲", prot, ig.get("亚洲")), ("蛋白×IgAN欧洲", prot, ig.get("欧洲")),
             ("蛋白×IGP29", prot, gly.get("IGP29")), ("IGP29×IgAN合并", gly.get("IGP29"), ig.get("合并")),
             ("IGP29×IgAN亚洲", gly.get("IGP29"), ig.get("亚洲")), ("IGP54×IgAN合并", gly.get("IGP54"), ig.get("合并"))]
    log("")
    log("B. 共定位先验敏感性（单因果变异ABF；表中为 PP3/PP4）")
    grid = [(1e-6, 0.15), (5e-6, 0.15), (1e-5, 0.15), (5e-5, 0.15), (1e-5, 0.10), (1e-5, 0.20)]
    log("  比较                " + "  ".join(f"p12={p:.0e},W={w}" for p, w in grid))
    for name, d1, d2 in comps:
        if not d1 or not d2:
            continue
        pr = pairs(d1, d2)
        if pr is None or len(pr[0]) < 50:
            continue
        _q, b1, s1, b2, s2 = pr
        cells = []
        for p12, w in grid:
            pp = coloc(b1, s1, b2, s2, p12=p12, w=w)
            cells.append(f"{pp[3]:.2f}/{pp[4]:.2f}")
            _rows.append(["prior", name, f"{p12:g}", f"{w}", f"{pp[3]:.4g}", f"{pp[4]:.4g}", len(b1), "", ""])
        log(f"  {name:<16} n={len(b1):<5}" + "        ".join(cells))

    # ---------- C. 条件共定位
    log("")
    log("C. 条件共定位：两边都先对一个信号做校正（z' = (z - r·z_k)/sqrt(1-r2)，剔除与该信号 r2>0.8 的位点）")

    def condition(d, kpos, bed):
        """d: {pos:(ea,oa,b,se)}，返回以bed A1为效应等位基因、对kpos校正后的 {pos:(A1,A2,b',se')}"""
        kk = key_of(kpos)
        if kk not in bed.idx or kpos not in d:
            return None
        ka1, ka2 = bed.idx[kk][1], bed.idx[kk][2]
        o = d[kpos]
        sk = m.sgn(o[0], o[1], ka1, ka2)
        if not sk:
            return None
        zk = sk * o[2] / o[3]
        gk = bed.geno(kk)[0]
        out = {}
        for q, o in d.items():
            kq = key_of(q)
            if q == kpos or kq not in bed.idx:
                continue
            a1, a2 = bed.idx[kq][1], bed.idx[kq][2]
            s = m.sgn(o[0], o[1], a1, a2)
            if not s:
                continue
            gq = bed.geno(kq)[0]
            ok = ~np.isnan(gq) & ~np.isnan(gk)
            if ok.sum() < 50 or gq[ok].std() == 0:
                continue
            r = float(np.corrcoef(gq[ok], gk[ok])[0, 1])
            if r * r > 0.8:
                continue
            z = s * o[2] / o[3]
            zc = (z - r * zk) / math.sqrt(1 - r * r)
            sec = o[3] / math.sqrt(1 - r * r)
            out[q] = (a1, a2, zc * sec, sec)
        return out
    plans = [("蛋白|S2 × IgAN欧洲|S2", prot, ig.get("欧洲"), S2, "EUR"), ("蛋白|S2 × IgAN合并|S2", prot, ig.get("合并"), S2, "EUR"),
             ("蛋白|S2 × IgAN亚洲|S2", prot, ig.get("亚洲"), S2, "EAS"),
             ("IGP29|S1 × IgAN欧洲|S1", gly.get("IGP29"), ig.get("欧洲"), S1, "EUR"), ("IGP29|S1 × IgAN合并|S1", gly.get("IGP29"), ig.get("合并"), S1, "EUR"),
             ("IGP29|S1 × IgAN亚洲|S1", gly.get("IGP29"), ig.get("亚洲"), S1, "EAS"),
             ("蛋白|S2 × IGP29|S2", prot, gly.get("IGP29"), S2, "EUR")]
    cache = {}
    for name, d1, d2, kpos, pop in plans:
        if not d1 or not d2 or pop not in beds:
            continue
        c1 = cache.setdefault((id(d1), kpos, pop), condition(d1, kpos, beds[pop]))
        c2 = cache.setdefault((id(d2), kpos, pop), condition(d2, kpos, beds[pop]))
        if not c1 or not c2:
            log(f"  {name}：校正信号在数据或参考面板中缺失，跳过")
            continue
        pr = pairs(c1, c2)
        if pr is None or len(pr[0]) < 50:
            log(f"  {name}：共有位点不足")
            continue
        q, b1, s1, b2, s2 = pr
        pp = coloc(b1, s1, b2, s2)
        z1, z2 = b1 / s1, b2 / s2
        i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
        log(f"  {name:<24} {pop} LD  位点 {len(q):,}  PP3={pp[3]:.2f} PP4={pp[4]:.2f}  "
            f"左侧最强 chr3:{q[i1]:,}（z={z1[i1]:+.1f}，右侧z={z2[i1]:+.1f}）；右侧最强 chr3:{q[i2]:,}（z={z2[i2]:+.1f}，左侧z={z1[i2]:+.1f}）")
        _rows.append(["conditional_coloc", name, pop, len(q), f"{pp[3]:.4g}", f"{pp[4]:.4g}", f"3:{q[i1]}", f"3:{q[i2]}", ""])


def main():
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    os.makedirs(os.path.dirname(RES), exist_ok=True)
    try:
        run()
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    with open(RES, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["analysis", "a", "b", "c", "d", "e", "f", "g", "h"])
        w.writerows(_rows)
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
