#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
122_olink_step1_egfrcys.py  —  第一步在Olink平台上的复制：eGFRcys → Olink蛋白（2026-10-04）

目的：稿件的新内容是"肾功能下降因果性地升高endotrophin"。Cho等2025已用CKDGen eGFR对
UKB Olink蛋白做过MR，但其附表拿不到，无法确认COL6A3是否在其中。这里直接用本项目已下载的
UKB-PPP全基因组汇总数据（raw\\protein\\ukbppp\\*.tar，37个蛋白）自己做一遍：
  暴露  Pan-UKB胱抑素C工具变量（464个）与CKDGen eGFRcys工具变量（159个），与主分析同一批
  结局  每个Olink蛋白（全基因组，trans位点为主）
  方法  与98步完全相同（import 98步的harmonise / ivw / weighted_median / mr_egger）
Olink的COL6A3抗体针对C端（endotrophin所在片段），所以这是对"肾功能→endotrophin"的跨平台复制。

位点对接：UKB-PPP文件的ID列是 chr:pos:a0:a1:imp:v1（GRCh37坐标），GENPOS是GRCh38。
工具变量表的POS是GRCh37，优先用ID里的坐标对；对不上的再用 logs\\_instruments_pos37_v4.tsv
里的GRCh38坐标对GENPOS。等位基因在harmonise里核对。

注意：Pan-UKB暴露与UKB-PPP蛋白同出UK Biobank（样本重叠），估计会偏向观察性关联；
CKDGen工具变量的重叠小一些。两套都报，稿件里如实说明。

COL6A3排第一个跑，结果每跑完一个蛋白就落盘，中途关掉重跑会跳过已完成的。
每个tar约550MB，预计每个蛋白1-3分钟，37个蛋白共约1小时；只想看COL6A3可以跑完第一个就关。

输出
  results\\olink_step1_egfrcys.tsv
  logs\\olink_step1_egfrcys_report.txt
"""

import csv
import glob
import gzip
import importlib.util
import io
import os
import sys
import tarfile
import time
from collections import Counter

import numpy as np

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
TAR_DIR = os.path.join(PROJ, "raw", "protein", "ukbppp")
OUT = os.path.join(RES, "olink_step1_egfrcys.tsv")
REPORT = os.path.join(LOG_DIR, "olink_step1_egfrcys_report.txt")
SETS = [("PanUKB", os.path.join(EXT, "panukb_cystatin_instruments.tsv")),
        ("CKDGen", os.path.join(EXT, "egfrcys_instruments.tsv"))]
POS38 = os.path.join(LOG_DIR, "_instruments_pos37_v4.tsv")
FIRST = ["col6a3", "clmp", "hs6st1", "apof", "fn1", "il6r", "gsn", "ptk7", "plg", "apoh"]
COLS = ["Gene", "instrument_set", "nSNP", "beta_IVW", "se_IVW", "p_IVW", "Q_p", "I2",
        "beta_WM", "se_WM", "p_WM", "beta_Egger", "se_Egger", "p_Egger", "egger_intercept_p",
        "beta_SomaScan", "p_SomaScan", "SomaScan_SeqIds", "n_matched_by_id37", "n_matched_by_pos38"]
_out = []


def log(m=""):
    print(m, flush=True)
    _out.append(str(m))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_sets():
    """返回 instr[set][rsid]=dict, key37{(chr,pos)->rsid}, key38{(chr,pos)->rsid}"""
    instr, key37, key38 = {}, {}, {}
    p38 = {}
    if os.path.exists(POS38):
        with open(POS38, "r", encoding="utf-8") as f:
            f.readline()
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) >= 5 and v[4].isdigit():
                    p38[v[2]] = (v[3].replace("chr", ""), v[4])
    for name, path in SETS:
        d = {}
        if not os.path.exists(path):
            log(f"!! 缺 {path}")
            instr[name] = d
            continue
        with open(path, "r", encoding="utf-8") as f:
            cols = f.readline().rstrip("\n").split("\t")
            for line in f:
                v = dict(zip(cols, line.rstrip("\n").split("\t")))
                d[v["SNP"]] = dict(EA=v["EA"].upper(), OA=v["OA"].upper(), BETA=float(v["BETA"]),
                                   SE=float(v["SE"]), EAF=float(v["EAF"]))
                key37[(str(v["CHR"]).replace("chr", ""), str(int(float(v["POS"]))))] = v["SNP"]
                if v["SNP"] in p38:
                    key38[p38[v["SNP"]]] = v["SNP"]
        instr[name] = d
        log(f"{name}: {len(d)} 个工具变量")
    log(f"位置键：GRCh37 {len(key37)} 个，GRCh38 {len(key38)} 个")
    return instr, key37, key38


def scan_tar(path, key37, key38):
    """返回 rows=[(rsid, ea, oa, beta, se)], 计数"""
    rows, seen, cnt = [], set(), Counter()
    with tarfile.open(path, "r") as tf:
        for member in tf:
            if not member.isfile():
                continue
            fobj = tf.extractfile(member)
            raw = gzip.GzipFile(fileobj=fobj) if member.name.endswith(".gz") else fobj
            stream = io.TextIOWrapper(io.BufferedReader(raw, 1 << 20), encoding="utf-8", errors="replace")
            header = stream.readline().split()
            low = [h.lower() for h in header]
            try:
                i_c, i_p, i_id = low.index("chrom"), low.index("genpos"), low.index("id")
                i_a0, i_a1 = low.index("allele0"), low.index("allele1")
                i_b, i_se = low.index("beta"), low.index("se")
            except ValueError:
                cnt["bad_header"] += 1
                if cnt["bad_header"] == 1:
                    log(f"    !! 表头对不上：{member.name}: {header[:14]}")
                continue
            mx = max(i_c, i_p, i_id, i_a0, i_a1, i_b, i_se)
            for line in stream:
                v = line.split(None, mx + 1)
                if len(v) <= mx:
                    continue
                idp = v[i_id].split(":", 2)
                rs = key37.get((idp[0].replace("chr", ""), idp[1])) if len(idp) >= 2 else None
                how = "id37"
                if rs is None:
                    rs = key38.get((v[i_c].replace("chr", ""), v[i_p]))
                    how = "pos38"
                if rs is None:
                    continue
                a0, a1 = v[i_a0].upper(), v[i_a1].upper()
                if len(a0) != 1 or len(a1) != 1:
                    continue
                try:
                    b, se = float(v[i_b]), float(v[i_se])
                except ValueError:
                    continue
                if (rs, a0, a1) in seen:
                    continue
                seen.add((rs, a0, a1))
                rows.append((rs, a1, a0, b, se))       # ALLELE1是效应等位基因
                cnt[how] += 1
    return rows, cnt


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    m = _load("_mod98", "98_mr_step1_panukb.py")
    instr, key37, key38 = load_sets()
    soma = {}
    p = os.path.join(RES, "mr_step1_panukb_to_protein.tsv")
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                soma.setdefault(r["Gene"].upper(), []).append(r)
    done = {}
    if os.path.exists(OUT):
        with open(OUT, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                done[(r["Gene"], r["instrument_set"])] = r
    tars = sorted(glob.glob(os.path.join(TAR_DIR, "*.tar")))
    if os.environ.get("MR_LIMIT"):
        tars = tars[:int(os.environ["MR_LIMIT"])]

    def rank(t):
        g = os.path.basename(t).split("_")[0].lower()
        return (FIRST.index(g) if g in FIRST else len(FIRST), g)
    tars.sort(key=rank)
    log(f"共 {len(tars)} 个Olink蛋白；已完成 {len({k[0] for k in done})} 个")
    results = list(done.values())
    for t in tars:
        base = os.path.basename(t)
        gene = base.split("_")[0].upper()
        if gene == "MICB":
            gene = "MICB_MICA"
        if all((gene, s) in done for s, _ in SETS):
            continue
        t0 = time.time()
        log("")
        log(f"--- {gene}  ({base})")
        try:
            rows, cnt = scan_tar(t, key37, key38)
        except Exception as e:  # noqa: BLE001
            log(f"    !! 读tar失败：{type(e).__name__}: {e}")
            continue
        log(f"    命中 {len(rows)} 行（按ID里的GRCh37坐标 {cnt['id37']}，按GENPOS {cnt['pos38']}），"
            f"用时 {time.time() - t0:.0f}s")
        ss = soma.get(gene.split("_")[0], [])
        for sname, _path in SETS:
            counters = Counter()
            bx, bxse, by, byse = m.harmonise(rows, instr[sname], counters)
            if len(bx) < 10:
                log(f"    {sname}: 可用工具变量只有 {len(bx)} 个，跳过（{dict(counters)}）")
                continue
            b, se, pv, _q, qp, i2 = m.ivw(bx, by, byse)
            wb, wse, wp = m.weighted_median(bx, bxse, by, byse)
            eg = m.mr_egger(bx, by, byse)
            eb, ese, ep, eip = (eg[0], eg[1], eg[2], eg[5]) if len(eg) == 6 else (np.nan,) * 4
            r = dict(Gene=gene, instrument_set=sname, nSNP=len(bx), beta_IVW=f"{b:.6g}", se_IVW=f"{se:.6g}",
                     p_IVW=f"{pv:.6g}", Q_p=f"{qp:.6g}", I2=f"{i2:.4g}", beta_WM=f"{wb:.6g}",
                     se_WM=f"{wse:.6g}", p_WM=f"{wp:.6g}", beta_Egger=f"{eb:.6g}", se_Egger=f"{ese:.6g}",
                     p_Egger=f"{ep:.6g}", egger_intercept_p=f"{eip:.6g}",
                     beta_SomaScan=";".join(x["beta_IVW"] for x in ss) if sname == "PanUKB" else "",
                     p_SomaScan=";".join(x["p_IVW"] for x in ss) if sname == "PanUKB" else "",
                     SomaScan_SeqIds=";".join(x["SeqId"] for x in ss),
                     n_matched_by_id37=cnt["id37"], n_matched_by_pos38=cnt["pos38"])
            results.append(r)
            log(f"    {sname}: n={len(bx)}  IVW beta={b:.4f} (SE {se:.4f}) P={pv:.3g}  "
                f"加权中位数 {wb:.4f} P={wp:.3g}  Egger截距P={eip:.3g}  "
                f"[协调：回文剔除{counters['palindromic']}，等位不符{counters['allele_mismatch']}]")
            if sname == "PanUKB" and ss:
                log("      SomaScan同基因：" + "；".join(
                    f"{x['SeqId']} beta={float(x['beta_IVW']):.3f} P={float(x['p_IVW']):.2g}" for x in ss))
        with open(OUT, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(results)
        with open(REPORT, "w", encoding="utf-8") as f:
            f.write("\n".join(_out) + "\n")
    log("")
    log("beta方向与主分析一致：每SD更高eGFRcys（Pan-UKB）或每单位log eGFRcys（CKDGen）对应的蛋白SD变化；"
        "负值=肾功能下降时蛋白升高")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    print(f"\n结果 {OUT}\n报告 {REPORT}")


if __name__ == "__main__":
    main()
