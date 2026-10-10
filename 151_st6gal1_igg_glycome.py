# -*- coding: utf-8 -*-
"""151步：ST6GAL1位点 × IgG N-糖组（Klaric 2020，77个性状，N约8,090）
对每个糖链性状：
  1) 蛋白最强pQTL处的Wald（每SD循环ST6GAL1对应的糖链性状变化）
  2) 蛋白 × 糖链 共定位
  3) 糖链 × IgA肾病（欧洲/亚洲/合并）共定位
第一次运行要把3 GB压缩包完整读一遍（约20-40分钟），窗口数据会缓存，之后再跑只要1分钟。
只读，不联网。
"""
import csv
import glob
import gzip
import io
import math
import os
import re
import sys
import tarfile
import time
import traceback
import importlib.util

import numpy as np

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "st6gal1_igg_glycome_report.txt")
GD = os.path.join(PROJ, "raw", "outcome", "igg_glycome_klaric2020")
GTAR = os.path.join(GD, "IgG_glycans_GWAS.tar.gz")
GMAP = os.path.join(GD, "map_r37.txt")
CACHE = os.path.join(EXT, "igg_glycome_st6gal1_window.tsv")
GENE, CHROM, LEAD_POS, LEAD_RS = "ST6GAL1", "3", 186739677, "rs4686837"
WIN = 500_000
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


def window_rsids():
    d = {}
    lo, hi = LEAD_POS - WIN, LEAD_POS + WIN
    with open(GMAP, "r", encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            v = line.split()
            if len(v) >= 3 and v[1] == CHROM:
                try:
                    p = int(v[2])
                except ValueError:
                    continue
                if lo <= p <= hi:
                    d[v[0]] = p
    return d


def read_glycans(rs2pos):
    """返回 {trait: {pos37: (ea, oa, beta, se)}}；有缓存就读缓存"""
    g = {}
    if os.path.exists(CACHE):
        with open(CACHE, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                g.setdefault(r["trait"], {})[int(r["pos37"])] = (r["ea"], r["oa"], float(r["beta"]), float(r["se"]))
        log(f"  读缓存 {CACHE}：{len(g)} 个性状")
        return g
    want = {k.encode(): v for k, v in rs2pos.items()}
    t0 = time.time()
    with tarfile.open(GTAR, "r|gz") as tf:
        for m in tf:
            if not m.isfile():
                continue
            mm = re.search(r"(IGP\d+)", m.name)
            if not mm:
                continue
            trait = mm.group(1)
            raw = gzip.GzipFile(fileobj=tf.extractfile(m))
            st = io.BufferedReader(raw, 1 << 22)
            head = st.readline().decode().lower().split()
            try:
                ie, io_, ib, ise = (head.index(x) for x in ("effect_allele", "other_allele", "effect", "se"))
            except ValueError:
                log(f"  !! 表头对不上 {m.name}: {head}")
                continue
            d = {}
            for line in st:
                rs = line[:line.find(b"\t")]
                pos = want.get(rs)
                if pos is None:
                    continue
                v = line.decode().split()
                try:
                    b, se = float(v[ib]), float(v[ise])
                except (ValueError, IndexError):
                    continue
                ea, oa = v[ie].upper(), v[io_].upper()
                if se > 0 and len(ea) == 1 and len(oa) == 1:
                    d[pos] = (ea, oa, b, se)
            g[trait] = d
            log(f"  {trait:<6} 窗口内 {len(d):,} 个位点   已用 {time.time() - t0:.0f}s")
    os.makedirs(EXT, exist_ok=True)
    with open(CACHE, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", "pos37", "ea", "oa", "beta", "se"])
        for t, d in g.items():
            for p, (ea, oa, b, se) in sorted(d.items()):
                w.writerow([t, p, ea, oa, b, se])
    return g


def pair(m, d1, d2, flip12):
    """d1: {pos:(...)}，返回对齐后的数组。flip12(pos)->(b1,s1,b2,s2)或None"""
    pos, b1, s1, b2, s2 = [], [], [], [], []
    for q in d1:
        if q not in d2:
            continue
        r = flip12(q)
        if r is None:
            continue
        pos.append(q); b1.append(r[0]); s1.append(r[1]); b2.append(r[2]); s2.append(r[3])
    return np.array(pos), np.array(b1), np.array(s1), np.array(b2), np.array(s2)


def run():
    m = load("_m147", "147_olink_igan_coloc.py")
    for p in (GTAR, GMAP):
        if not os.path.exists(p) and not os.path.exists(CACHE):
            log(f"找不到 {p}")
            return
    log("A. 糖组数据：取ST6GAL1领头pQTL ±500 kb")
    rs2pos = window_rsids() if os.path.exists(GMAP) else {}
    log(f"  map里窗口内位点 {len(rs2pos):,} 个；领头 {LEAD_RS} " + ("在" if LEAD_RS in rs2pos else "不在") + "糖组数据的位点表里")
    gly = read_glycans(rs2pos)
    if not gly:
        log("  没读到糖组数据")
        return

    log("")
    log("B. Olink ST6GAL1 与 IgA肾病")
    tars = [t for t in glob.glob(os.path.join(m.TAR_DIR, "*.tar")) if os.path.basename(t).upper().startswith(GENE + "_")]
    if not tars:
        log("  找不到ST6GAL1的Olink tar")
        return
    ol = m.read_olink(tars[0], CHROM, LEAD_POS - WIN, LEAD_POS + WIN)
    log(f"  Olink窗口内 {len(ol):,} 个位点")
    ig = {}
    for name, path in m.IGAN:
        if os.path.exists(path):
            d, _n = m.scan_igan(path, {CHROM: [(LEAD_POS - WIN, LEAD_POS + WIN)]})
            ig[name] = {k[1]: v for k, v in d.items()}
            log(f"  IgA肾病 {name}：窗口内 {len(ig[name]):,} 个位点")

    log("")
    log("C. 逐个糖链性状（Wald=每SD循环ST6GAL1对应的性状变化，性状已做秩正态化；按P排序）")
    res = []
    for trait in sorted(gly, key=lambda t: int(t[3:])):
        gd = gly[trait]

        def f1(q, gd=gd):
            a0, a1, b, se = ol[q]
            ea, oa, gb, gse = gd[q]
            s = m.sgn(a1, a0, ea, oa)
            return None if s == 0 else (b, se, s * gb, gse)
        pos, b1, s1, b2, s2 = pair(m, ol, gd, f1)
        if len(pos) < 50:
            log(f"  {trait}：与Olink共有位点只有 {len(pos)} 个，跳过")
            continue
        pp = m.coloc(b1, s1, b2, s2)
        z1, z2 = b1 / s1, b2 / s2
        i1, i2 = int(np.argmax(np.abs(z1))), int(np.argmax(np.abs(z2)))
        wald, wse = b2[i1] / b1[i1], abs(s2[i1] / b1[i1])
        row = dict(trait=trait, n=len(pos), wald=wald, wse=wse, p=pz(wald / wse), pp3=pp[3], pp4=pp[4],
                   top_pqtl=int(pos[i1]), minp=pz(z2[i2]), top_gly=int(pos[i2]), top_gly_pz=z1[i2])
        for name, idat in ig.items():
            def f2(q, gd=gd, idat=idat):
                ea, oa, gb, gse = gd[q]
                o = idat[q]
                s = m.sgn(ea, oa, o[0], o[1])
                return None if s == 0 else (gb, gse, s * o[2], o[3])
            p2, c1, e1, c2, e2 = pair(m, gd, idat, f2)
            row["igan_" + name] = float(m.coloc(c1, e1, c2, e2)[4]) if len(p2) >= 50 else float("nan")
        res.append(row)
    res.sort(key=lambda r: r["p"])
    names = list(ig)
    log("  性状    位点  Wald (se)          P        蛋白×糖链 PP3/PP4   糖链窗口最小P（该位点蛋白z）   糖链×IgA肾病PP4（" + "/".join(names) + "）")
    for r in res:
        igs = " / ".join(f"{r['igan_' + n]:.2f}" for n in names)
        log(f"  {r['trait']:<6} {r['n']:>5}  {r['wald']:+.3f} ({r['wse']:.3f})  {r['p']:.1e}   {r['pp3']:.2f}/{r['pp4']:.2f}   "
            f"{r['minp']:.1e}（z={r['top_gly_pz']:+.1f}）   {igs}")
        _rows.append([r["trait"], r["n"], f"{r['wald']:.5g}", f"{r['wse']:.5g}", f"{r['p']:.4g}", f"{r['pp3']:.4g}", f"{r['pp4']:.4g}",
                      f"{CHROM}:{r['top_pqtl']}", f"{r['minp']:.4g}", f"{CHROM}:{r['top_gly']}", f"{r['top_gly_pz']:.2f}"]
                     + [f"{r['igan_' + n]:.4g}" for n in names])
    log("")
    thr = 0.05 / max(len(res), 1)
    hit = [r for r in res if r["p"] < thr]
    both = [r for r in hit if r["pp4"] >= 0.8]
    log(f"D. 汇总：{len(res)} 个性状；Wald P<{thr:.1e}（Bonferroni）的 {len(hit)} 个，其中蛋白×糖链PP4>=0.8的 {len(both)} 个")
    if res:
        log(f"   蛋白最强pQTL（与糖组共有位点中）：chr{CHROM}:{res[0]['top_pqtl']:,}；领头 {LEAD_RS} 在 chr{CHROM}:{LEAD_POS:,}")
    for r in both:
        tri = [n for n in names if r["igan_" + n] >= 0.8]
        log(f"   {r['trait']}: Wald {r['wald']:+.3f}，PP4 {r['pp4']:.2f}；与IgA肾病共定位（PP4>=0.8）的数据集：{','.join(tri) if tri else '无'}")
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "st6gal1_igg_glycome.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["trait", "n_snps", "wald", "se", "p", "PP3", "PP4", "top_pQTL", "glycan_min_p", "top_glycan", "top_glycan_pQTL_z"]
                   + ["PP4_glycan_IgAN_" + n for n in names])
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
