#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
116_phewas_finemap_eqtl.py  —  对标Nature Genetics范本（Yoshiji 2025）补三块公开数据证据
（2026-10-04）

  A. 领头pQTL的全表型扫描（PheWAS）：rs2020854（APOF）、rs11677932（COL6A3 C端）
     从 FinnGen、Biobank Japan PheWeb、UKB-TOPMed PheWeb、UKB-SAIGE PheWeb 的公开API
     各取一次JSON（每个几百KB），原样存盘。等位基因顺序两种都试。
  B. GTEx v8：两个基因的各组织中位表达；领头pQTL在相关组织里是不是该基因的eQTL
     （dyneqtl接口，任意变异-基因-组织都能查，返回效应量NES和P）。
  C. 精细定位（本地计算，不联网）：APOF、COL6A3在Aragam CAD和蛋白两侧，
     把其余独立信号条件掉之后对主信号做ABF 95%可信集，报告领头pQTL的后验概率(PIP)
     以及两侧可信集的重叠。复用104步的 Ref / joint_condition / independent_signals。
  D. APOF主信号 × 115步已缓存的13个性状：只留APOF主pQTL信号（条件掉其余9个）后，
     在±100kb窄窗里重算共定位——回答115步"全区域共定位不可解读"留下的问题。
     只读 raw\\panukb_slices\\ 里的缓存，不重新下载。

A、B两部分只负责下载并原样保存JSON（接口字段各家不同，解析放到之后做），
任何一个接口失败都不影响其余部分。

输出
  raw\\phewas\\*.json            各接口原始返回
  results\\finemap_credible_sets.tsv
  results\\apof_main_signal_trait_coloc.tsv
  logs\\phewas_finemap_eqtl_report.txt
"""

import importlib.util
import json
import math
import os
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request

import numpy as np

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.join(PROJ, "extracted")
RES = os.path.join(PROJ, "results")
LOG_DIR = os.path.join(PROJ, "logs")
PHE_DIR = os.path.join(PROJ, "raw", "phewas")
SLICE_DIR = os.path.join(PROJ, "raw", "panukb_slices")
REPORT = os.path.join(LOG_DIR, "phewas_finemap_eqtl_report.txt")
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-data-download"
SKIP_NET = os.environ.get("MR_SKIP_NET") == "1"

# rsid, 基因, GRCh38 (chr,pos), GRCh37 (chr,pos), 两个等位基因
VARIANTS = [("rs2020854", "APOF", ("12", 56349583), ("12", 56743367), ("T", "C")),
            ("rs11677932", "COL6A3", ("2", 237315312), ("2", 238223955), ("G", "A"))]
PHEWEB = [("finngen_r10", "https://r10.finngen.fi/api/variant/{c}-{p}-{a}-{b}", 38),
          ("finngen_r12", "https://r12.finngen.fi/api/variant/{c}-{p}-{a}-{b}", 38),
          ("bbj", "https://pheweb.jp/api/variant/{c}-{p}-{a}-{b}", 37),
          ("ukb_topmed", "https://pheweb.org/UKB-TOPMed/api/variant/{c}-{p}-{a}-{b}", 38),
          ("ukb_saige", "https://pheweb.org/UKB-SAIGE/api/variant/{c}-{p}-{a}-{b}", 37)]
GTEX = "https://gtexportal.org/api/v2/"
GTEX_TISSUES = {"APOF": ["Liver", "Kidney_Cortex", "Whole_Blood", "Adipose_Subcutaneous",
                         "Artery_Coronary", "Artery_Aorta"],
                "COL6A3": ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial",
                           "Cells_Cultured_fibroblasts", "Adipose_Subcutaneous",
                           "Adipose_Visceral_Omentum", "Kidney_Cortex", "Liver", "Whole_Blood"]}
GENES = [("12370_30", "APOF"), ("11196_31", "COL6A3")]
ARAGAM = os.path.join(EXT, "coloc_panukb", "CAD.tsv")
_out = []


def log(m=""):
    print(m, flush=True)
    _out.append(str(m))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def fetch_json(url, dest, tries=3):
    """下载并原样保存；返回解析后的对象或None。已存在就直接读。"""
    if os.path.exists(dest) and os.path.getsize(dest) > 2:
        try:
            with open(dest, "r", encoding="utf-8") as f:
                return json.load(f)
        except ValueError:
            pass
    if SKIP_NET:
        return None
    last = ""
    for i in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                raw = r.read()
            obj = json.loads(raw.decode("utf-8", "replace"))
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "wb") as f:
                f.write(raw)
            return obj
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
            if e.code in (400, 403, 404, 422):
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError, ValueError) as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(3 * i)
    log(f"      失败（{last}）：{url}")
    return None


def find_records(obj):
    """在任意嵌套JSON里找"最像表型结果"的那个列表（元素是dict且带p值字段）。"""
    best = []

    def walk(o):
        nonlocal best
        if isinstance(o, list):
            if o and isinstance(o[0], dict) and any(
                    k in o[0] for k in ("pval", "pvalue", "pValue", "p")):
                if len(o) > len(best):
                    best = o
            for x in o[:50]:
                walk(x)
        elif isinstance(o, dict):
            for v in o.values():
                walk(v)
    walk(obj)
    return best


# ------------------------------------------------------------------ A
def part_a():
    log("A. 领头pQTL全表型扫描（原始JSON存到 raw\\phewas\\）")
    for rs, gene, hg38, hg37, (a1, a2) in VARIANTS:
        log(f"  {rs}（{gene}）")
        for name, tmpl, build in PHEWEB:
            c, p = hg38 if build == 38 else hg37
            got = None
            for x, y in ((a1, a2), (a2, a1)):
                url = tmpl.format(c=c, p=p, a=x, b=y)
                dest = os.path.join(PHE_DIR, f"{name}_{rs}_{x}_{y}.json")
                obj = fetch_json(url, dest)
                recs = find_records(obj) if obj is not None else []
                if recs:
                    got = (x, y, recs)
                    break
            if not got:
                log(f"    {name:<12} 没取到结果")
                continue
            x, y, recs = got

            def pv(r):
                for k in ("pval", "pvalue", "pValue", "p"):
                    try:
                        return float(r[k])
                    except (KeyError, TypeError, ValueError):
                        continue
                return 1.0
            recs = sorted(recs, key=pv)
            log(f"    {name:<12} ref-alt按 {x}-{y} 取到 {len(recs)} 个表型；P最小的8个：")
            for r in recs[:8]:
                nm = r.get("phenostring") or r.get("phenotype") or r.get("phenocode") or r.get("pheno") or "?"
                b = r.get("beta", r.get("effect", ""))
                log(f"        P={pv(r):.2g}  beta={b}  {str(nm)[:70]}")


# ------------------------------------------------------------------ B
def part_b():
    log("")
    log("B. GTEx v8（原始JSON存到 raw\\phewas\\gtex_*.json）")
    for rs, gene, hg38, _hg37, (a1, a2) in VARIANTS:
        g = fetch_json(GTEX + "reference/gene?" + urllib.parse.urlencode({"geneId": gene}),
                       os.path.join(PHE_DIR, f"gtex_gene_{gene}.json"))
        gid = None
        for r in find_any(g, "gencodeId"):
            gid = r
            break
        if not gid:
            log(f"  {gene}：没取到gencodeId，跳过")
            continue
        log(f"  {gene}：gencodeId={gid}")
        ex = fetch_json(GTEX + "expression/medianGeneExpression?" +
                        urllib.parse.urlencode({"gencodeId": gid}),
                        os.path.join(PHE_DIR, f"gtex_expr_{gene}.json"))
        rows = [r for r in (ex or {}).get("data", []) if isinstance(r, dict)]
        rows.sort(key=lambda r: -float(r.get("median", 0) or 0))
        if rows:
            log("    中位表达最高的6个组织：" + "；".join(
                f"{r.get('tissueSiteDetailId')}={float(r.get('median', 0)):.1f}" for r in rows[:6]))
            kid = [r for r in rows if "Kidney" in str(r.get("tissueSiteDetailId"))]
            for r in kid:
                log(f"    {r.get('tissueSiteDetailId')}={float(r.get('median', 0)):.2f} TPM")
        c, p = hg38
        for tis in GTEX_TISSUES[gene]:
            done = False
            for x, y in ((a1, a2), (a2, a1)):
                vid = f"chr{c}_{p}_{x}_{y}_b38"
                d = fetch_json(GTEX + "association/dyneqtl?" + urllib.parse.urlencode(
                    {"gencodeId": gid, "variantId": vid, "tissueSiteDetailId": tis}),
                    os.path.join(PHE_DIR, f"gtex_eqtl_{gene}_{rs}_{tis}_{x}_{y}.json"))
                if isinstance(d, dict) and d.get("pValue") is not None:
                    log(f"    {rs} × {gene} 在 {tis:<28} NES={d.get('nes')}  P={d.get('pValue')}"
                        f"  （{vid}，NES对应alt等位基因）")
                    done = True
                    break
            if not done:
                log(f"    {rs} × {gene} 在 {tis:<28} 没取到")


def find_any(obj, key):
    if isinstance(obj, dict):
        if key in obj and isinstance(obj[key], str):
            yield obj[key]
        for v in obj.values():
            yield from find_any(v, key)
    elif isinstance(obj, list):
        for x in obj:
            yield from find_any(x, key)


# ------------------------------------------------------------------ C / D 公用
def abf_pip(z, se, w=0.15 ** 2):
    v = se ** 2
    r = w / (w + v)
    labf = 0.5 * (np.log(1 - r) + r * z ** 2)
    labf = labf - labf.max()
    p = np.exp(labf)
    return p / p.sum()


def credible(pip, names, level=0.95):
    order = np.argsort(-pip)
    cum, cs = 0.0, []
    for i in order:
        cs.append(names[i])
        cum += pip[i]
        if cum >= level:
            break
    return cs


def setup_locus(m, ref, prot, outc):
    rows = m.harmonise(prot, outc)
    rsids = [r[0] for r in rows]
    R, keep = ref.ld_matrix(rsids, {rs: prot[rs][0] for rs in rsids})
    pos = {rs: i for i, rs in enumerate(keep)}
    sel = [r for r in rows if r[0] in pos]
    order = [pos[r[0]] for r in sel]
    R = R[np.ix_(order, order)]
    names = [r[0] for r in sel]
    b1 = np.array([r[1] for r in sel]); s1 = np.array([r[2] for r in sel])
    b2 = np.array([r[3] for r in sel]); s2 = np.array([r[4] for r in sel])
    bp = np.array([prot[rs][4] for rs in names], float)
    return names, R, b1, s1, b2, s2, bp


def part_c(m, ref, prots):
    log("")
    log("C. 精细定位：条件掉其余信号后，主信号的ABF 95%可信集（欧洲LD，Aragam CAD）")
    if not os.path.exists(ARAGAM):
        log(f"  找不到 {ARAGAM}，跳过")
        return
    outc = m.read_outcome(ARAGAM)
    lead_of = {g: rs for rs, g, *_ in VARIANTS}
    rows_out = []
    for seq, gene in GENES:
        prot = prots[(seq, gene)]
        names, R, b1, s1, b2, s2, bp = setup_locus(m, ref, prot, outc)
        if len(names) < m.MIN_SNP:
            log(f"  {gene}：位点不足，跳过")
            continue
        ix = {rs: i for i, rs in enumerate(names)}
        lead = lead_of[gene]
        if lead not in ix:
            log(f"  {gene}：领头位点{lead}不在可用位点里，跳过")
            continue
        z1, z2 = b1 / s1, b2 / s2
        leads1 = m.independent_signals(z1, R)
        others = [i for i in leads1 if R[i, ix[lead]] ** 2 < 0.5 and i != ix[lead]]
        zc1, ok1 = m.joint_condition(z1, R, others)
        use = ok1.copy(); use[ix[lead]] = True
        win = np.abs(bp - bp[ix[lead]]) <= 250_000
        u = use & win
        nm = [names[i] for i in np.where(u)[0]]
        pip1 = abf_pip(zc1[u], s1[u])
        pip2 = abf_pip(z2[u], s2[u])
        cs1, cs2 = credible(pip1, nm), credible(pip2, nm)
        li = nm.index(lead)
        inter = sorted(set(cs1) & set(cs2))
        top1 = nm[int(np.argmax(pip1))]; top2 = nm[int(np.argmax(pip2))]
        shared = float(np.sum(pip1 * pip2))
        log(f"  {gene}（{seq}）：±250kb内 {len(nm)} 个位点，蛋白侧条件掉 {len(others)} 个其他信号")
        log(f"    蛋白侧：{lead} PIP={pip1[li]:.3f}；95%可信集 {len(cs1)} 个位点，PIP最高 {top1}")
        log(f"    CAD侧 ：{lead} PIP={pip2[li]:.3f}；95%可信集 {len(cs2)} 个位点，PIP最高 {top2}"
            f"（与{lead}的r²={R[ix[top2], ix[lead]]**2:.3f}）")
        log(f"    两侧可信集交集 {len(inter)} 个：{inter[:12]}{'…' if len(inter) > 12 else ''}")
        log(f"    Σ(PIP蛋白×PIP_CAD)={shared:.3f}（两侧后验逐位点乘积之和，描述性指标，越接近1说明两侧指向同一位点）")
        rows_out.append(dict(Gene=gene, SeqId=seq, lead=lead, n_snps=len(nm),
                             n_protein_signals_conditioned=len(others),
                             pip_protein_lead=pip1[li], pip_cad_lead=pip2[li],
                             cs95_protein_n=len(cs1), cs95_cad_n=len(cs2),
                             cs95_overlap_n=len(inter), top_protein=top1, top_cad=top2,
                             r2_topcad_lead=R[ix[top2], ix[lead]] ** 2,
                             sum_pip_product=shared, cs95_overlap=",".join(inter)))
    cols = ["Gene", "SeqId", "lead", "n_snps", "n_protein_signals_conditioned",
            "pip_protein_lead", "pip_cad_lead", "cs95_protein_n", "cs95_cad_n",
            "cs95_overlap_n", "top_protein", "top_cad", "r2_topcad_lead", "sum_pip_product",
            "cs95_overlap"]
    with open(os.path.join(RES, "finemap_credible_sets.tsv"), "w", encoding="utf-8",
              newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows_out:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], (float, np.floating))
                              else str(r[c]) for c in cols) + "\n")


def part_d(m, ref, prots, m111):
    log("")
    log("D. APOF主pQTL信号（条件掉其余蛋白信号）× 115步缓存的性状：±100kb窄窗共定位")
    if not os.path.isdir(SLICE_DIR):
        log("  没有 raw\\panukb_slices（115步没跑过），跳过")
        return
    pos37 = m111.load_pos37(m111.PANUKB_INSTR)
    seq, gene, lead = "12370_30", "APOF", "rs2020854"
    prot = prots[(seq, gene)]
    rows_out = []
    files = sorted(f for f in os.listdir(SLICE_DIR) if f.endswith(".tsv") and f"_{gene}_chr" in f)
    for fn in files:
        trait = fn.split(f"_{gene}_chr")[0]
        outc = {}
        with open(os.path.join(SLICE_DIR, fn), "r", encoding="utf-8") as f:
            hd = f.readline().rstrip("\n").split("\t")
            hi = {c: i for i, c in enumerate(hd)}
            if any(c not in hi for c in ("chr", "pos", "ref", "alt", "beta_EUR", "se_EUR")):
                continue
            lc = hi.get("low_confidence_EUR")
            for line in f:
                v = line.rstrip("\n").split("\t")
                if len(v) < len(hd) or (lc is not None and v[lc].lower() == "true"):
                    continue
                try:
                    b, se = float(v[hi["beta_EUR"]]), float(v[hi["se_EUR"]])
                except ValueError:
                    continue
                if not (se > 0) or not math.isfinite(b):
                    continue
                for rs in pos37.get(f"{v[hi['chr']]}:{v[hi['pos']]}", []):
                    if rs in prot:
                        outc[rs] = (v[hi["alt"]].upper(), v[hi["ref"]].upper(), b, se)
        names, R, b1, s1, b2, s2, bp = setup_locus(m, ref, prot, outc)
        if len(names) < m.MIN_SNP or lead not in names:
            log(f"  {trait:<11} 可用位点不足或缺领头位点，跳过")
            continue
        ix = {rs: i for i, rs in enumerate(names)}
        z1, z2 = b1 / s1, b2 / s2
        leads1 = m.independent_signals(z1, R)
        others = [i for i in leads1 if R[i, ix[lead]] ** 2 < 0.5 and i != ix[lead]]
        zc1, ok1 = m.joint_condition(z1, R, others)
        use = ok1.copy(); use[ix[lead]] = True
        u = use & (np.abs(bp - bp[ix[lead]]) <= m.WINDOW_BP)
        if u.sum() < m.MIN_SNP_WINDOW:
            log(f"  {trait:<11} 窄窗位点不足（{int(u.sum())}），跳过")
            continue
        pp = m.coloc(zc1[u] * s1[u], s1[u], b2[u], s2[u])
        zl = z2[ix[lead]]
        wmin = pz(float(np.abs(z2[u]).max()))
        log(f"  {trait:<11} {lead}上该性状P={pz(zl):.2g}；窄窗内性状最小P={wmin:.2g}；"
            f"PP1={pp[1]:.3f} PP3={pp[3]:.3f} PP4={pp[4]:.3f}（{int(u.sum())}个位点）")
        rows_out.append(dict(trait=trait, n=int(u.sum()), p_trait_at_lead=pz(zl),
                             window_min_p=wmin, PP0=pp[0], PP1=pp[1], PP2=pp[2], PP3=pp[3],
                             PP4=pp[4]))
    cols = ["trait", "n", "p_trait_at_lead", "window_min_p", "PP0", "PP1", "PP2", "PP3", "PP4"]
    with open(os.path.join(RES, "apof_main_signal_trait_coloc.tsv"), "w", encoding="utf-8",
              newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows_out:
            f.write("\t".join(f"{r[c]:.6g}" if isinstance(r[c], (float, np.floating))
                              else str(r[c]) for c in cols) + "\n")


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    for d in (RES, LOG_DIR, PHE_DIR):
        os.makedirs(d, exist_ok=True)
    try:
        for fn, label in ((part_a, "A"), (part_b, "B")):
            try:
                fn()
            except Exception as e:  # noqa: BLE001
                log(f"!! {label}部分异常：{type(e).__name__}: {e}")
        try:
            m = _load("_mod104", "104_coloc_sensitivity_panukb.py")
            m111 = _load("_mod111", "111_independent_outcome_replication.py")
            if os.environ.get("MR_OLD_PROJ"):          # 仅测试用
                old = os.environ["MR_OLD_PROJ"]
                m.OLD_CIS_DIR = os.path.join(old, "extracted", "cis")
                m.LD_PREFIX = os.path.join(old, "raw", "reference", "ld", "EUR")
            prots = {(s, g): m.read_cis(s, g) for s, g in GENES}
            wanted = set()
            for d in prots.values():
                wanted |= set(d)
            ref = m.Ref(m.LD_PREFIX, wanted)
            for fn, label in ((lambda: part_c(m, ref, prots), "C"),
                              (lambda: part_d(m, ref, prots, m111), "D")):
                try:
                    fn()
                except Exception as e:  # noqa: BLE001
                    log(f"!! {label}部分异常：{type(e).__name__}: {e}")
                    for ln in traceback.format_exc().splitlines():
                        log("    " + ln)
        except Exception as e:  # noqa: BLE001
            log(f"!! C/D准备阶段异常：{type(e).__name__}: {e}")
            for ln in traceback.format_exc().splitlines():
                log("    " + ln)
    finally:
        with open(REPORT, "w", encoding="utf-8") as f:
            f.write("\n".join(_out) + "\n")
        print(f"\n报告已写入 {REPORT}")


if __name__ == "__main__":
    main()
