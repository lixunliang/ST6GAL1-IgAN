#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
146_download_olink_igan_candidates.py  —  下载IgA肾病候选蛋白的Olink全量数据（2026-10-07）

145步用Olink领头位点筛出了一批IgA肾病候选，但只有领头位点做不了共定位。这里从Synapse取这些蛋白的
UKB-PPP欧洲发现集全量数据（每个约550 MB），下一步做共定位、次级信号检验和位点图。
复用131步的下载函数（读 logs\\ukbppp_file_listing.tsv，不重新列目录；已有的文件跳过）。
需要环境变量 SYNAPSE_AUTH_TOKEN（由bat提示粘贴，不落盘）。

候选（25个）：TNFSF13(APRIL) TNFRSF13B(TACI) TNFRSF13C(BAFF-R) CFH CFHR5 CFB C3 ST6GAL1 ITGAM PYDC1 FCAR
              FCGR2A FCGR2B FCRL2 FCRL3 PADI4 EIF2AK3 NFKB1 IL7R CD28 RNASET2 SAT2 TNFRSF6B CXCL6 PIGR
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GENES = ["TNFSF13", "TNFRSF13B", "TNFRSF13C", "CFH", "CFHR5", "CFB", "C3", "ST6GAL1", "ITGAM", "PYDC1", "FCAR",
         "FCGR2A", "FCGR2B", "FCRL2", "FCRL3", "PADI4", "EIF2AK3", "NFKB1", "IL7R", "CD28", "RNASET2", "SAT2",
         "TNFRSF6B", "CXCL6", "PIGR"]

if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    spec = importlib.util.spec_from_file_location("_m131", os.path.join(HERE, "131_download_kidney_line_data.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.GENES = GENES
    m.do_olink()
