#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
98_mr_step1_panukb.py  —  第一步 MR（Pan-UKB版本）：胱抑素C定义的肾功能 → 4,907 个血浆蛋白

跟 30_mr_step1_egfrcys.py 统计方法完全一致（IVW随机效应/加权中位数/MR-Egger/
Cochran's Q/I²/BH-FDR），唯一区别是暴露端换成 95 步产出的 Pan-UKB 胱抑素C
工具变量（464个，R²≈14.06%，是CKDGen那159个的~70倍），配上96步b扫出来的
对应deCODE蛋白效应。两套结果不是互相替代，是两个独立数据源的triangulation——
CAD和LIFESPAN这两个总效应只有Pan-UKB版本才显著（97步已验证），这一步就是要看
这两条通路背后到底是哪些蛋白在起中介作用。

方向：panukb_cystatin_instruments.tsv 在95步已经把原始"胱抑素C浓度"方向翻转
成跟 egfrcys_instruments.tsv 一致的"beta越大=肾功能越好"，所以这里的
beta_IVW 解读规则跟30步完全一样：beta>0=蛋白随'肾功能变好'而升高，
beta<0=蛋白随'肾功能下降'而升高（更贴近"下游损伤标志物"的直觉）。

输入
  extracted\\panukb_cystatin_instruments.tsv       464 个 Pan-UKB 胱抑素C 工具变量（95步产出）
  extracted\\step1_panukb_cystatin\\shard_*.tsv    96步b 扫出来的各蛋白在这些位点上的效应（deCODE）

输出
  results\\mr_step1_panukb_to_protein.tsv          每个蛋白一行的完整结果
  results\\mr_step1_panukb_harmonisation.txt       协调过程的计数（写 Methods 用）
"""

import argparse
import glob
import os
import sys
import time

import numpy as np

try:
    from scipy import stats
except ImportError:  # noqa: BLE001
    sys.exit("缺少 scipy，请在 Anaconda Prompt 里跑：conda install scipy")

PROJ = r"F:\孟德尔"
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
INSTR = os.path.join(EXT, "panukb_cystatin_instruments.tsv")
SHARDS = os.path.join(EXT, "step1_panukb_cystatin", "shard_*.tsv")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROMIC = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
N_BOOT = 500
RNG = np.random.default_rng(20260928)


def log(m):
    print(m, flush=True)


# ------------------------------------------------------------------ 读入
def load_instruments():
    d = {}
    with open(INSTR, "r", encoding="utf-8") as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = dict(zip(cols, line.rstrip("\n").split("\t")))
            d[v["SNP"]] = dict(EA=v["EA"].upper(), OA=v["OA"].upper(),
                               BETA=float(v["BETA"]), SE=float(v["SE"]),
                               EAF=float(v["EAF"]))
    return d


def load_shards(instr):
    """返回 {seqid: (gene, [(rsid, ea, oa, beta, se), ...])}"""
    prot = {}
    n_rows = 0
    for path in sorted(glob.glob(SHARDS)):
        with open(path, "r", encoding="utf-8") as f:
            cols = f.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(cols)}
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) < len(cols):
                    continue
                n_rows += 1
                rs = v[ix["rsids"]]
                if rs not in instr:
                    if "," in rs:
                        hit = [x for x in rs.split(",") if x in instr]
                        if not hit:
                            continue
                        rs = hit[0]
                    else:
                        continue
                seq = v[ix["SeqId"]]
                g = prot.setdefault(seq, (v[ix["Gene"]], []))
                try:
                    g[1].append((rs, v[ix["effectAllele"]].upper(),
                                 v[ix["otherAllele"]].upper(),
                                 float(v[ix["Beta"]]), float(v[ix["SE"]])))
                except ValueError:
                    continue
    return prot, n_rows


# ------------------------------------------------------------------ 协调
def harmonise(rows, instr, counters):
    bx, bxse, by, byse = [], [], [], []
    seen = set()
    for rs, dea, doa, dbeta, dse in rows:
        if rs in seen:
            continue
        e = instr.get(rs)
        if e is None:
            continue
        if (e["EA"], e["OA"]) in PALINDROMIC:
            counters["palindromic"] += 1
            seen.add(rs)
            continue
        if dea not in COMP or doa not in COMP:
            counters["indel_or_other"] += 1
            continue
        sign = 0
        if (dea, doa) == (e["EA"], e["OA"]):
            sign = 1
        elif (dea, doa) == (e["OA"], e["EA"]):
            sign = -1
        else:
            cdea, cdoa = COMP[dea], COMP[doa]
            if (cdea, cdoa) == (e["EA"], e["OA"]):
                sign = 1
                counters["strand_flipped"] += 1
            elif (cdea, cdoa) == (e["OA"], e["EA"]):
                sign = -1
                counters["strand_flipped"] += 1
        if sign == 0:
            counters["allele_mismatch"] += 1
            continue
        if dse <= 0 or e["SE"] <= 0:
            continue
        seen.add(rs)
        bx.append(e["BETA"])
        bxse.append(e["SE"])
        by.append(sign * dbeta)
        byse.append(dse)
        counters["kept"] += 1
    return (np.array(bx), np.array(bxse), np.array(by), np.array(byse))


# ------------------------------------------------------------------ MR
def ivw(bx, by, byse):
    w = 1.0 / byse**2
    denom = np.sum(w * bx**2)
    beta = np.sum(w * bx * by) / denom
    se_fe = np.sqrt(1.0 / denom)
    resid = by - beta * bx
    q = float(np.sum(w * resid**2))
    df = len(bx) - 1
    phi = max(1.0, q / df) if df > 0 else 1.0
    se = se_fe * np.sqrt(phi)
    p = 2 * stats.norm.sf(abs(beta / se))
    q_p = float(stats.chi2.sf(q, df)) if df > 0 else float("nan")
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 and df > 0 else 0.0
    return beta, se, p, q, q_p, i2


def weighted_median_estimate(ratio, w):
    order = np.argsort(ratio)
    r, ww = ratio[order], w[order]
    cw = np.cumsum(ww) - 0.5 * ww
    cw /= np.sum(ww)
    k = np.searchsorted(cw, 0.5)
    if k == 0:
        return r[0]
    if k >= len(r):
        return r[-1]
    below, above = cw[k - 1], cw[k]
    if above == below:
        return r[k]
    return r[k - 1] + (r[k] - r[k - 1]) * (0.5 - below) / (above - below)


def weighted_median(bx, bxse, by, byse):
    ratio = by / bx
    vr = byse**2 / bx**2 + (by**2 * bxse**2) / bx**4
    w = 1.0 / vr
    est = weighted_median_estimate(ratio, w)
    bxb = RNG.normal(bx, bxse, size=(N_BOOT, len(bx)))
    byb = RNG.normal(by, byse, size=(N_BOOT, len(by)))
    with np.errstate(divide="ignore", invalid="ignore"):
        rb = byb / bxb
    boots = np.empty(N_BOOT)
    for i in range(N_BOOT):
        ok = np.isfinite(rb[i])
        boots[i] = weighted_median_estimate(rb[i][ok], w[ok]) if ok.any() else np.nan
    se = float(np.nanstd(boots, ddof=1))
    p = 2 * stats.norm.sf(abs(est / se)) if se > 0 else float("nan")
    return est, se, p


def mr_egger(bx, by, byse):
    n = len(bx)
    if n < 4:
        return (np.nan,) * 5
    w = 1.0 / byse**2
    X = np.column_stack([np.ones(n), bx])
    W = np.diag(w)
    try:
        xtwx_inv = np.linalg.inv(X.T @ W @ X)
    except np.linalg.LinAlgError:
        return (np.nan,) * 5
    coef = xtwx_inv @ (X.T @ W @ by)
    resid = by - X @ coef
    df = n - 2
    sigma2 = float((resid * w * resid).sum() / df)
    cov = xtwx_inv * max(sigma2, 1.0)   # 与 TwoSampleMR 一致：过离散不缩小 SE
    se = np.sqrt(np.diag(cov))
    t_slope = coef[1] / se[1]
    t_int = coef[0] / se[0]
    p_slope = 2 * stats.t.sf(abs(t_slope), df)
    p_int = 2 * stats.t.sf(abs(t_int), df)
    return coef[1], se[1], p_slope, coef[0], se[0], p_int


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    pp = p[ok]
    n = len(pp)
    order = np.argsort(pp)
    ranked = pp[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


# ------------------------------------------------------------------ 主流程
def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-snp", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(RES, exist_ok=True)
    t0 = time.time()

    if not os.path.exists(INSTR):
        sys.exit(f"缺少工具变量文件 {INSTR}，请先跑 95_build_panukb_cystatin_instruments.py")
    instr = load_instruments()
    log(f"Pan-UKB 胱抑素C 工具变量 {len(instr):,} 个")

    shard_paths = sorted(glob.glob(SHARDS))
    if not shard_paths:
        sys.exit(f"缺少扫描分片 {SHARDS}，请先跑 96步b（96_scan_decode_panukb_cystatin.py --workers 5）")
    prot, n_rows = load_shards(instr)
    log(f"读入 {n_rows:,} 行、{len(prot):,} 个蛋白")

    counters = dict(palindromic=0, strand_flipped=0, allele_mismatch=0,
                    indel_or_other=0, kept=0)
    out = []
    for i, (seq, (gene, rows)) in enumerate(sorted(prot.items()), 1):
        bx, bxse, by, byse = harmonise(rows, instr, counters)
        n = len(bx)
        if n < args.min_snp:
            continue
        b, se, p, q, qp, i2 = ivw(bx, by, byse)
        wb, wse, wp = weighted_median(bx, bxse, by, byse)
        eb, ese, ep, ei, eise, eip = mr_egger(bx, by, byse)
        out.append(dict(SeqId=seq, Gene=gene, nSNP=n,
                        beta_IVW=b, se_IVW=se, p_IVW=p,
                        Q=q, Q_p=qp, I2=i2,
                        beta_WM=wb, se_WM=wse, p_WM=wp,
                        beta_Egger=eb, se_Egger=ese, p_Egger=ep,
                        egger_intercept=ei, egger_intercept_se=eise,
                        egger_intercept_p=eip))
        if i % 500 == 0:
            log(f"  已算 {i:,}/{len(prot):,} 个蛋白 ...")

    ps = [r["p_IVW"] for r in out]
    qs = bh_fdr(ps)
    for r, qv in zip(out, qs):
        r["p_FDR"] = qv

    cols = ["SeqId", "Gene", "nSNP", "beta_IVW", "se_IVW", "p_IVW", "p_FDR",
            "Q", "Q_p", "I2", "beta_WM", "se_WM", "p_WM",
            "beta_Egger", "se_Egger", "p_Egger",
            "egger_intercept", "egger_intercept_se", "egger_intercept_p"]
    res_path = os.path.join(RES, "mr_step1_panukb_to_protein.tsv")
    out.sort(key=lambda r: r["p_IVW"])
    with open(res_path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(
                (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]))
                for c in cols) + "\n")

    sig = [r for r in out if r["p_FDR"] < 0.05]
    robust = [r for r in sig
              if np.isfinite(r["p_WM"]) and r["p_WM"] < 0.05
              and np.sign(r["beta_WM"]) == np.sign(r["beta_IVW"])
              and (not np.isfinite(r["egger_intercept_p"])
                   or r["egger_intercept_p"] >= 0.05)]
    up = sum(1 for r in sig if r["beta_IVW"] > 0)

    summary = [
        "第一步 MR（Pan-UKB版：胱抑素C定义的肾功能 → 血浆蛋白）协调与结果汇总",
        f"分析蛋白数：{len(out):,} / {len(prot):,}（工具变量数 < {args.min_snp} 的已剔除）",
        f"每个蛋白平均可用工具变量：{np.mean([r['nSNP'] for r in out]):.1f}",
        "",
        "方向提醒：beta>0 = 该蛋白随'肾功能变好'而升高；",
        "         beta<0 = 该蛋白随'肾功能下降'而升高（通常更符合'下游损伤标志物'的直觉）。",
        "（Pan-UKB原始方向是胱抑素C浓度，95步已翻转成跟eGFRcys一致的方向，可直接对照CKDGen版结果）",
        "",
        "等位基因协调计数（累计跨所有蛋白）：",
        f"  保留                {counters['kept']:,}",
        f"  回文位点剔除        {counters['palindromic']:,}",
        f"  按互补链对齐        {counters['strand_flipped']:,}",
        f"  等位基因对不上剔除  {counters['allele_mismatch']:,}",
        f"  非 SNP/其他         {counters['indel_or_other']:,}",
        "",
        f"FDR<0.05 的蛋白：{len(sig):,}（随肾功能变好而升高 {up:,}，降低 {len(sig)-up:,}）",
        f"其中同时满足「加权中位数同向且 P<0.05」且「MR-Egger 截距 P>=0.05」"
        f"的稳健蛋白：{len(robust):,}",
        "",
        f"结果表：{res_path}",
        f"用时 {(time.time()-t0)/60:.1f} 分钟",
        "",
        "FDR<0.05 且方法间稳健、随'肾功能下降'而升高的前 25 个（按 P 排序）：",
    ]
    down_robust = [r for r in robust if r["beta_IVW"] < 0]
    for r in down_robust[:25]:
        summary.append(
            f"  {r['Gene']:<12} {r['SeqId']:<10} beta={r['beta_IVW']:+.3f} "
            f"P={r['p_IVW']:.2e} FDR={r['p_FDR']:.2e} nSNP={r['nSNP']} I2={r['I2']:.0f}%")

    txt = "\n".join(summary)
    with open(os.path.join(RES, "mr_step1_panukb_harmonisation.txt"), "w",
              encoding="utf-8") as f:
        f.write(txt + "\n")
    log("")
    log(txt)


if __name__ == "__main__":
    main()
