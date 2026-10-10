#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
131_download_kidney_line_data.py  —  肾功能线补充数据下载（2026-10-04）

  winkler  CKDGen糖尿病分层eGFR（Winkler 2022, Commun Biol）：一个文件里同时有
           糖尿病人群和非糖尿病人群的效应（Effect.dm / Effect.nodm）
           -> raw\\outcome\\ckdgen_winkler2022_dmstrat\\
  olink    UKB-PPP（Olink）欧洲发现集里肾功能线候选蛋白的汇总数据，按
           logs\\ukbppp_file_listing.tsv 里已有的清单直接取，不再重新列目录
           -> raw\\protein\\ukbppp\\
           三角验证的19个蛋白里Olink平台有9个：GSTA1 INHBB CTSS MST1 MICB BTD HEXIM1 AK1 F2
           （MANBA、POR、TCEA2、SFTA2、STAT3、DNAJB4、BMP7、CNPY3、FAIM、EFNA3不在Olink上）

用法  python 131_download_kidney_line_data.py winkler
      <小环境python> 131_download_kidney_line_data.py olink   （需要环境变量SYNAPSE_AUTH_TOKEN）
"""
import csv
import importlib.util
import os
import sys

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PROJ, "logs")
WINKLER_URL = "https://ckdgen.imbi.uni-freiburg.de/files/Winkler2022/meta_egfr_dmstrat_stage1plus2.txt.gz"
WINKLER_DIR = os.path.join(PROJ, "raw", "outcome", "ckdgen_winkler2022_dmstrat")
OLINK_DIR = os.path.join(PROJ, "raw", "protein", "ukbppp")
LISTING = os.path.join(LOG_DIR, "ukbppp_file_listing.tsv")
GENES = ["GSTA1", "INHBB", "CTSS", "MST1", "MICB", "BTD", "HEXIM1", "AK1", "F2"]
ANCESTRY = "European (discovery)"


def do_winkler():
    spec = importlib.util.spec_from_file_location(
        "_mod110", os.path.join(HERE, "110_download_independent_outcomes.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    report = os.path.join(LOG_DIR, "download_winkler2022_report.txt")
    m._out.clear()
    os.makedirs(WINKLER_DIR, exist_ok=True)
    m.log(f"CKDGen 糖尿病分层eGFR  {WINKLER_URL}")
    try:
        pk = m.peek(WINKLER_URL)
    except Exception as e:  # noqa: BLE001
        pk = None
        m.log(f"  !! 探表头异常：{type(e).__name__}: {e}")
    if pk is None:
        m.log("  !! 探表头失败（地址不通或被拦），没有下载")
    else:
        is_gz, header, samples, total = pk
        m.log(f"  期望大小：{total:,}字节" if total else "  期望大小未知")
        m.log(f"  表头：{header.split()}")
        for s in samples[:3]:
            m.log(f"  样例：{s[:400]}")
        ok = m.download(WINKLER_URL, os.path.join(WINKLER_DIR, os.path.basename(WINKLER_URL)), total, is_gz)
        m.log("  下载完成" if ok else "  !! 下载未成功（.part已保留，重跑会续传）")
    with open(report, "w", encoding="utf-8") as f:
        f.write("\n".join(m._out) + "\n")


def do_olink():
    out = []

    def log(s=""):
        print(s, flush=True)
        out.append(s)
    report = os.path.join(LOG_DIR, "download_olink_kidney_line_report.txt")
    try:
        if not os.path.exists(LISTING):
            log(f"!! 找不到 {LISTING}，请先跑100步a")
            return
        hits = []
        with open(LISTING, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["type"] != "FileEntity" or ANCESTRY not in r["path"]:
                    continue
                for g in GENES:
                    if r["name"].upper().startswith(g + "_"):
                        hits.append((r["synapse_id"], r["name"], g))
        log(f"清单里匹配到 {len(hits)} 个文件（应为 {len(GENES)} 个）")
        for sid, name, g in hits:
            log(f"  {g:<8} {sid:<13} {name}")
        token = os.environ.get("SYNAPSE_AUTH_TOKEN", "").strip()
        if not token:
            log("!! 没读到 SYNAPSE_AUTH_TOKEN，没有下载")
            return
        try:
            import synapseclient
        except ImportError:
            log("!! 这个python里没装synapseclient，请用bat运行")
            return
        os.makedirs(OLINK_DIR, exist_ok=True)
        syn = synapseclient.Synapse()
        syn.login(authToken=token, silent=True)
        log("Synapse 登录成功")
        for i, (sid, name, g) in enumerate(hits, 1):
            target = os.path.join(OLINK_DIR, name)
            if os.path.exists(target) and os.path.getsize(target) > 0:
                log(f"  [{i}/{len(hits)}] {name} 已存在（{os.path.getsize(target):,}字节），跳过")
                continue
            try:
                syn.get(sid, downloadLocation=OLINK_DIR, ifcollision="keep.local")
                sz = os.path.getsize(target) if os.path.exists(target) else -1
                log(f"  [{i}/{len(hits)}] {g} {name} 完成（{sz:,}字节）")
            except Exception as e:  # noqa: BLE001
                log(f"  [{i}/{len(hits)}] {name} 失败：{type(e).__name__}: {e}")
        log("全部结束")
    finally:
        os.makedirs(LOG_DIR, exist_ok=True)
        with open(report, "w", encoding="utf-8") as f:
            f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "winkler":
        do_winkler()
    elif mode == "olink":
        do_olink()
    else:
        print("用法：131_download_kidney_line_data.py winkler|olink")
