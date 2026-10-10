# -*- coding: utf-8 -*-
"""159步：外周血单细胞数据 GSE285335（健康9、早期IgA肾病11、晚期6）里B细胞的ST6GAL1
做法（不依赖单细胞软件，规则全部写明）：
  1) 逐个样本下载10x矩阵（每个40-85 MB，共约1.6 GB，存 raw\\geo\\GSE285335\\，可断点重跑）
  2) 用标志基因划细胞：B细胞 = CD79A>=1 且 (MS4A1>=1 或 CD19>=1) 且 CD3E=0 且 CD14=0 且 LYZ<=1
     记忆B = B细胞且CD27>=1；初始B = B细胞且CD27=0且(IGHD>=1或TCL1A>=1)
     浆母/浆细胞 = JCHAIN>=3 且 MZB1>=1 且 MS4A1=0 且 CD3E=0；T细胞 = CD3E>=1 且 CD79A=0 且 CD14=0；单核 = (CD14>=1或LYZ>=3) 且 CD3E=0 且 CD79A=0
  3) 每个样本、每类细胞做"拟bulk"：该类细胞里目标基因UMI之和 / 该类细胞总UMI × 1e6（CPM），以及表达该基因的细胞比例
  4) 以样本为单位比较三组（Mann-Whitney，正态近似）
联网，约30-60分钟（主要是下载）。
"""
import csv
import gzip
import math
import os
import sys
import time
import traceback
import urllib.request

import numpy as np

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
D = os.path.join(PROJ, "raw", "geo", "GSE285335")
RES = os.path.join(PROJ, "results")
REPORT = os.path.join(PROJ, "logs", "scrna_pbmc_st6gal1_report.txt")
BASE = os.environ.get("GEO_BASE", "https://ftp.ncbi.nlm.nih.gov/geo/samples")
SAMPLES = [(8700986 + i, f"L{i + 1}") for i in range(6)] + [(8700992 + i, f"E{i + 1}") for i in range(11)] + \
          [(8701003 + i, f"H{i + 1}") for i in range(9)]
MARK = ["CD79A", "MS4A1", "CD19", "CD3E", "CD14", "LYZ", "CD27", "IGHD", "TCL1A", "JCHAIN", "MZB1"]
TARGET = ["ST6GAL1", "CD22", "B4GALT1", "ST3GAL1", "C1GALT1", "C1GALT1C1", "ST6GALNAC2", "TNFRSF13B", "TNFRSF13C", "ACTB"]
TYPES = ["B细胞", "记忆B", "初始B", "浆母细胞", "T细胞", "单核"]
GRP = {"H": "健康", "E": "早期IgAN", "L": "晚期IgAN"}
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        return True
    last = ""
    for i in range(4):
        try:
            if url.startswith("file:"):
                import shutil
                shutil.copy(url[5:], dest)
                return True
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=300) as r, open(dest + ".part", "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            os.replace(dest + ".part", dest)
            return True
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
            time.sleep(5 * (i + 1))
    log(f"    下载失败：{last[:160]}  {url}")
    return False


def read_mtx(path):
    """返回 (行号数组, 列号数组, 数值数组, 基因数, 细胞数)，行列从0开始"""
    try:
        import pandas as pd
        with gzip.open(path, "rt") as f:
            skip = 0
            for line in f:
                skip += 1
                if not line.startswith("%"):
                    ng, nc, _nnz = (int(x) for x in line.split())
                    break
        df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=skip, dtype=np.float64, engine="c")
        return df[0].to_numpy().astype(np.int64) - 1, df[1].to_numpy().astype(np.int64) - 1, df[2].to_numpy(), ng, nc
    except ImportError:
        rr, cc, vv = [], [], []
        with gzip.open(path, "rt") as f:
            first = True
            for line in f:
                if line.startswith("%"):
                    continue
                a = line.split()
                if first:
                    ng, nc = int(a[0]), int(a[1])
                    first = False
                    continue
                rr.append(int(a[0]) - 1); cc.append(int(a[1]) - 1); vv.append(float(a[2]))
        return np.array(rr), np.array(cc), np.array(vv), ng, nc


def one(gsm, tag):
    pre = f"GSM{gsm}"
    stub = pre[:-3] + "nnn"
    fm = os.path.join(D, f"{pre}_{tag}_matrix.mtx.gz")
    ff = os.path.join(D, f"{pre}_features.tsv.gz")
    if not fetch(f"{BASE}/{stub}/{pre}/suppl/{pre}_{tag}_matrix.mtx.gz", fm):
        return None
    if not fetch(f"{BASE}/{stub}/{pre}/suppl/{pre}_features.tsv.gz", ff):
        return None
    sym = []
    with gzip.open(ff, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            v = line.rstrip("\n").split("\t")
            sym.append(v[1] if len(v) > 1 else v[0])
    r, c, v, ng, nc = read_mtx(fm)
    if ng != len(sym):
        log(f"    {tag}: 基因数对不上（矩阵 {ng}，features {len(sym)}）")
    tot = np.bincount(c, weights=v, minlength=nc)
    ngene = np.bincount(c, minlength=nc)
    idx = {}
    for i, s in enumerate(sym):
        idx.setdefault(s, i)
    g = {}
    for s in MARK + TARGET:
        x = np.zeros(nc)
        if s in idx:
            m = r == idx[s]
            np.add.at(x, c[m], v[m])
        g[s] = x
    ok = (tot >= 500) & (ngene >= 200)
    B = ok & (g["CD79A"] >= 1) & ((g["MS4A1"] >= 1) | (g["CD19"] >= 1)) & (g["CD3E"] == 0) & (g["CD14"] == 0) & (g["LYZ"] <= 1)
    sets = {"B细胞": B, "记忆B": B & (g["CD27"] >= 1), "初始B": B & (g["CD27"] == 0) & ((g["IGHD"] >= 1) | (g["TCL1A"] >= 1)),
            "浆母细胞": ok & (g["JCHAIN"] >= 3) & (g["MZB1"] >= 1) & (g["MS4A1"] == 0) & (g["CD3E"] == 0),
            "T细胞": ok & (g["CD3E"] >= 1) & (g["CD79A"] == 0) & (g["CD14"] == 0),
            "单核": ok & ((g["CD14"] >= 1) | (g["LYZ"] >= 3)) & (g["CD3E"] == 0) & (g["CD79A"] == 0)}
    rows = []
    for t, m in sets.items():
        n = int(m.sum())
        for s in TARGET:
            cpm = float(g[s][m].sum() / tot[m].sum() * 1e6) if n and tot[m].sum() > 0 else float("nan")
            frac = float((g[s][m] > 0).mean()) if n else float("nan")
            rows.append([tag, tag[0], t, n, int(ok.sum()), s, cpm, frac])
    log(f"  {tag:<4} 细胞 {nc:,}（过滤后 {int(ok.sum()):,}）  B {int(B.sum()):,}  记忆B {int(sets['记忆B'].sum()):,}  浆母 {int(sets['浆母细胞'].sum()):,}  T {int(sets['T细胞'].sum()):,}  单核 {int(sets['单核'].sum()):,}")
    return rows


def mwu(a, b):
    a, b = [x for x in a if not math.isnan(x)], [x for x in b if not math.isnan(x)]
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    allv = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    ranks, i = [0.0] * len(allv), 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        for k in range(i, j + 1):
            ranks[k] = (i + j) / 2 + 1
        i = j + 1
    ra = sum(rk for rk, (_x, w) in zip(ranks, allv) if w == 0)
    u = ra - len(a) * (len(a) + 1) / 2
    mu, sd = len(a) * len(b) / 2, math.sqrt(len(a) * len(b) * (len(a) + len(b) + 1) / 12)
    return math.erfc(abs(u - mu) / sd / math.sqrt(2)) if sd else float("nan")


def main():
    os.makedirs(D, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    cache = os.path.join(RES, "scrna_pbmc_st6gal1_pseudobulk.tsv")
    allrows, done = [], set()
    if os.path.exists(cache):
        with open(cache, "r", encoding="utf-8") as f:
            for r in csv.reader(f, delimiter="\t"):
                if r and r[0] != "sample":
                    allrows.append([r[0], r[1], r[2], int(r[3]), int(r[4]), r[5], float(r[6]), float(r[7])])
                    done.add(r[0])
    log("A. 逐个样本")
    for gsm, tag in SAMPLES:
        if tag in done:
            continue
        try:
            rows = one(gsm, tag)
        except Exception as e:  # noqa: BLE001
            log(f"  {tag}: 异常 {type(e).__name__}: {e}")
            for ln in traceback.format_exc().splitlines()[-4:]:
                log("      " + ln)
            rows = None
        if rows:
            allrows += rows
            with open(cache, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["sample", "group", "celltype", "n_cells", "n_cells_total", "gene", "cpm", "frac_expressing"])
                w.writerows(allrows)
    got = sorted(set(r[0] for r in allrows))
    log(f"  已完成样本 {len(got)}/{len(SAMPLES)}")
    log("")
    log("B. 各类细胞占比（占过滤后细胞的百分比，组内中位数）")
    for t in TYPES:
        txt = []
        for gk, gname in GRP.items():
            vals = [100 * r[3] / r[4] for r in allrows if r[2] == t and r[1] == gk and r[5] == TARGET[0] and r[4] > 0]
            txt.append(f"{gname} {np.median(vals):.2f}%（n={len(vals)}）" if vals else f"{gname} -")
        log(f"  {t:<6} " + "  ".join(txt))
    log("")
    log("C. 拟bulk表达（CPM，组内中位数；括号内为表达该基因的细胞比例的中位数）；P为Mann-Whitney，以样本为单位；细胞数<20的样本不计")
    for s in TARGET:
        log(f"  {s}")
        for t in TYPES:
            by, fr = {}, {}
            for r in allrows:
                if r[5] == s and r[2] == t and r[3] >= 20:
                    by.setdefault(r[1], []).append(r[6]); fr.setdefault(r[1], []).append(r[7])
            if not by:
                continue
            txt = "  ".join(f"{GRP[k]} {np.median(by[k]):.1f}（{100 * np.median(fr[k]):.0f}%，n={len(by[k])}）" for k in GRP if k in by)
            ig = by.get("E", []) + by.get("L", [])
            p1, p2, p3 = mwu(ig, by.get("H", [])), mwu(by.get("E", []), by.get("H", [])), mwu(by.get("L", []), by.get("H", []))
            log(f"    {t:<6} {txt}   P：全部IgAN对健康 {p1:.2g}；早期对健康 {p2:.2g}；晚期对健康 {p3:.2g}")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
