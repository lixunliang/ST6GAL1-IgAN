#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
143_download_finngen_kidney_diseases.py  —  下载FinnGen R12的肾脏疾病终点（2026-10-06）

新方向：蛋白 → 具体肾脏疾病（首选IgA肾病）。主数据是Kiryluk 2023的IgA肾病GWAS（10,146例，
需用户到哥伦比亚大学Kiryluk实验室网页手动下载）。FinnGen作独立复制和扩展病种：
  N14_IGA_NEPHROPATHY   IgA肾病               822例
  DM_NEPHROPATHY        糖尿病肾病          5,579例
  N14_GLOMER_NEPHRITIS  肾小球肾炎          3,905例
  N14_NEPHROTICSYND     肾病综合征          1,106例
  N14_GLOMERULAR        肾小球疾病（总）    7,426例
病例数来自142步读到的FinnGen R12终点清单。复用110步的探表头、断点续传、gzip校验。
需要联网，每个文件约700到800 MB。
输出  raw\\outcome\\finngen_r12_<终点小写>\\ ；logs\\download_finngen_kidney_diseases_report.txt
"""
import importlib.util
import os
import sys

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(PROJ, "raw", "outcome")
REPORT = os.path.join(PROJ, "logs", "download_finngen_kidney_diseases_report.txt")
FG = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{}.gz"
ENDPOINTS = ["N14_IGA_NEPHROPATHY", "DM_NEPHROPATHY", "N14_GLOMER_NEPHRITIS", "N14_NEPHROTICSYND", "N14_GLOMERULAR"]


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    spec = importlib.util.spec_from_file_location("_mod110", os.path.join(HERE, "110_download_independent_outcomes.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m._out.clear()
    for ep in ENDPOINTS:
        key = "finngen_r12_" + ep.lower()
        outdir = os.path.join(OUT_DIR, key)
        os.makedirs(outdir, exist_ok=True)
        pk, url = None, ""
        for url in (FG.format(ep), FG.format(ep).replace("/release/", "/")):
            m.log("=" * 60)
            m.log(f"{key}  {url}")
            try:
                pk = m.peek(url)
            except Exception as e:  # noqa: BLE001
                pk = None
                m.log(f"  !! 探表头异常：{type(e).__name__}: {e}")
            if pk is not None:
                break
        if pk is None:
            m.log("  !! 探表头失败，跳过")
            continue
        is_gz, header, samples, total = pk
        m.log(f"  期望大小：{total:,}字节" if total else "  期望大小未知")
        m.log(f"  表头：{header.split()}")
        if not ("beta" in header.lower() and "pval" in header.lower()):
            m.log("  !! 表头不像FinnGen汇总数据，不下载")
            continue
        if not m.download(url, os.path.join(outdir, os.path.basename(url)), total, is_gz):
            m.log("  !! 下载未成功（.part已保留，重跑会续传）")
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(m._out) + "\n")


if __name__ == "__main__":
    main()
