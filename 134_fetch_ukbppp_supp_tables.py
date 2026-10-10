#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
134_fetch_ukbppp_supp_tables.py  —  取UKB-PPP（Sun 2023, Nature）附表里的全蛋白组pQTL清单（2026-10-05）

目的：把肾功能线从"SomaScan筛选 + Olink验证9个"扩成"双平台全蛋白组"。Olink约2,900个蛋白的
全量汇总数据有1.6TB，不现实；论文附表里有每个蛋白的领头pQTL（含cis/trans、效应、标准误），
足够做全蛋白组的cis-MR筛选，筛出来的再按需下载单个蛋白的全量数据做共定位。

这一步只做下载和"看表头"：依次尝试Nature附件的几个候选地址，下载成功的xlsx逐个工作表打印
表名、行列数、前4行。我看过真实表头后再写筛选脚本（地址和表结构我没法事先核实）。

输出  raw\\reference\\ukbppp_sun2023\\ ；logs\\ukbppp_supp_tables_report.txt
"""
import os
import sys
import time
import traceback
import urllib.error
import urllib.request
import zipfile

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
OUT = os.path.join(PROJ, "raw", "reference", "ukbppp_sun2023")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "ukbppp_supp_tables_report.txt")
BASE = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/"
CAND = [f"41586_2023_6592_MOESM{i}_ESM.{ext}" for i in range(1, 9) for ext in ("xlsx", "zip")]
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/126 Safari/537.36"
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def fetch(name):
    dest = os.path.join(OUT, name)
    if os.path.exists(dest) and os.path.getsize(dest) > 10000:
        log(f"  {name} 已存在（{os.path.getsize(dest):,}字节）")
        return dest
    url = BASE + name
    for i in range(2):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r, open(dest + ".part", "wb") as f:
                n = 0
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
                    n += len(b)
                    if n % (20 << 20) < (1 << 20):
                        print(f"      {n / 1e6:.0f} MB", flush=True)
            os.replace(dest + ".part", dest)
            log(f"  {name} 下载完成（{os.path.getsize(dest):,}字节）")
            return dest
        except urllib.error.HTTPError as e:
            if e.code in (403, 404):
                return None
            log(f"  {name} HTTP {e.code}，重试")
        except Exception as e:  # noqa: BLE001
            log(f"  {name} {type(e).__name__}: {e}")
        time.sleep(3)
    return None


def peek_xlsx(path):
    try:
        import openpyxl
    except ImportError:
        log("    !! 这个python没装openpyxl，无法看表头（文件已下好；可运行 pip install openpyxl 后重跑）")
        return
    try:
        wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    except Exception as e:  # noqa: BLE001
        log(f"    !! 打不开：{type(e).__name__}: {e}")
        return
    for ws in wb.worksheets:
        log(f"    工作表 [{ws.title}]  约 {ws.max_row} 行 × {ws.max_column} 列")
        for k, row in enumerate(ws.iter_rows(values_only=True)):
            if k >= 6:
                break
            cells = ["" if c is None else str(c)[:28] for c in row[:40]]
            while cells and cells[-1] == "":
                cells.pop()
            if cells:
                log("      " + " | ".join(cells))
    wb.close()


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    try:
        got = []
        log("尝试下载 Sun 2023（Nature 622:329）附件：")
        for name in CAND:
            p = fetch(name)
            if p:
                got.append(p)
        # 手动放进来的文件也一并看
        for fn in sorted(os.listdir(OUT)):
            p = os.path.join(OUT, fn)
            if p not in got and fn.lower().endswith((".xlsx", ".zip")):
                got.append(p)
        if not got:
            log("!! 一个都没下到（地址不对或被拦）。请手动从文章页面 Supplementary Information 下载")
            log("   Supplementary Tables（xlsx），放进 raw\\reference\\ukbppp_sun2023\\ 后重跑本步。")
            log("   文章：https://www.nature.com/articles/s41586-023-06592-6")
        for p in got:
            log("")
            log(f"== {os.path.basename(p)}  {os.path.getsize(p):,}字节")
            if p.lower().endswith(".zip"):
                try:
                    with zipfile.ZipFile(p) as z:
                        for zi in z.infolist():
                            log(f"    {zi.filename}  {zi.file_size:,}")
                            if zi.filename.lower().endswith(".xlsx"):
                                z.extract(zi, OUT)
                                peek_xlsx(os.path.join(OUT, zi.filename))
                except Exception as e:  # noqa: BLE001
                    log(f"    !! zip打不开：{e}")
            else:
                peek_xlsx(p)
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")


if __name__ == "__main__":
    main()
