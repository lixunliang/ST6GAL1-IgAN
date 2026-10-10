#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
142_probe_kidney_disease_gwas.py  —  新方向探路：有哪些"肾脏疾病"GWAS可用（2026-10-06）

前面的结局都是肾功能指标（eGFR、尿素氮、蛋白尿）。要为国自然打基础，结局应换成具体的肾脏疾病
（首选IgA肾病，其次糖尿病肾病、膜性肾病等），这样筛出的蛋白才能接到机制研究上。
这一步只探路、不下大文件：
  1 从128步存下的FinnGen R12全表型结果里，列出所有肾脏相关终点及病例数（不联网）
  2 向GWAS Catalog查询几篇肾脏疾病GWAS的登记信息和是否有全量数据（联网，复用135步的函数）：
     Kiryluk 2023  IgA肾病（Nat Genet）                PMID 36797366
     Xie 2020      膜性肾病（Nat Commun）              PMID 32231244
     Salem 2019    糖尿病肾病（JASN，DNCRI）           PMID 31537649
     Sakaue 2021   日本人群多病种（含IgA肾病等）        PMID 34594039（只列肾脏相关条目）
PMID是凭记忆写的，报告会打印文章标题供核对。
输出 logs\\probe_kidney_disease_gwas_report.txt
"""
import glob
import importlib.util
import json
import os
import re
import sys
import traceback

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(PROJ, "logs", "probe_kidney_disease_gwas_report.txt")
KW = re.compile(r"kidney|renal|nephr|glomerul|iga |iga$|dialys|uraemi|uremi|proteinuri|albuminuri|"
                r"diabetic.*(kidney|nephro)|polycystic|lupus neph|tubulo", re.I)
PMIDS = [("Kiryluk 2023 IgA肾病", "36797366", None), ("Xie 2020 膜性肾病", "32231244", None),
         ("Salem 2019 糖尿病肾病", "31537649", None), ("Sakaue 2021 日本多病种", "34594039", KW)]


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    spec = importlib.util.spec_from_file_location("_m135", os.path.join(HERE, "135_probe_exposure_gwas.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    log = m.log
    try:
        log("1. FinnGen R12 肾脏相关终点（来自128步保存的全表型结果）")
        files = sorted(glob.glob(os.path.join(PROJ, "raw", "phewas", "finngen_r12_*.json")))
        rows = []
        if files:
            with open(files[0], "r", encoding="utf-8") as f:
                obj = json.load(f)

            def walk(o):
                if isinstance(o, list):
                    for x in o:
                        walk(x)
                elif isinstance(o, dict):
                    if "phenocode" in o and "phenostring" in o:
                        rows.append(o)
                    else:
                        for v in o.values():
                            walk(v)
            walk(obj)
        seen = set()
        hit = []
        for r in rows:
            if r["phenocode"] in seen:
                continue
            seen.add(r["phenocode"])
            if KW.search(str(r.get("phenostring", ""))) or str(r["phenocode"]).startswith("N14_"):
                hit.append(r)
        hit.sort(key=lambda r: -(r.get("n_case") or 0))
        log(f"  共 {len(seen)} 个终点，肾脏相关 {len(hit)} 个：")
        for r in hit:
            log(f"    {str(r['phenocode']):<38} 病例 {int(r.get('n_case') or 0):>7}  对照 {int(r.get('n_control') or 0):>7}  {str(r.get('phenostring'))[:70]}")
        if not files:
            log("  !! 没找到 raw\\phewas\\finngen_r12_*.json")
        log("")
        log("2. GWAS Catalog")
        for label, pmid, flt in PMIDS:
            log("=" * 70)
            log(f"{label}  PMID {pmid}")
            raw = m.get(m.API.format(pmid))
            if raw is None:
                continue
            try:
                studies = json.loads(raw.decode("utf-8", "replace")).get("_embedded", {}).get("studies", [])
            except ValueError:
                log("    返回的不是JSON")
                continue
            if studies:
                pi = studies[0].get("publicationInfo") or {}
                log(f"  文章：{(pi.get('author') or {}).get('fullname', '')} {pi.get('publicationDate', '')} "
                    f"{pi.get('publication', '')}｜{str(pi.get('title', ''))[:110]}")
            log(f"  登记研究 {len(studies)} 个")
            for s in studies:
                trait = (s.get("diseaseTrait") or {}).get("trait", "")
                if flt is not None and not flt.search(trait):
                    continue
                acc = s.get("accessionId", "")
                full = s.get("fullPvalueSet")
                log(f"  {acc:<14} 全量数据={'有' if full else '无'}  {trait[:55]:<55} | {str(s.get('initialSampleSize', ''))[:130]}")
                if full:
                    ents = m.list_dir(f"{m.FTP}{m.bucket(acc)}/{acc}/")
                    log(f"      目录：{ents if ents is not None else '列不出来'}")
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(m._out) + "\n")


if __name__ == "__main__":
    main()
