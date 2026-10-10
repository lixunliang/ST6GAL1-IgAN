# -*- coding: utf-8 -*-
"""158步：公共外周血表达数据里IgA肾病与对照的ST6GAL1
A. 芯片数据（series matrix，小文件）：GSE73953（PBMC，IgAN/膜性/对照）、GSE125818（PBMC）、GSE58539（单核细胞）
   同时列出B细胞标志基因（CD19、MS4A1、CD79A），用来判断差异是否只是B细胞比例不同
B. 单细胞数据 GSE285335（PBMC，健康/早期/晚期）：这一步只列出附件清单和大小，不下载
联网，约5分钟。文件存到 raw\\geo\\。
"""
import gzip
import io
import math
import os
import re
import sys
import time
import traceback
import urllib.request

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
GEO = os.path.join(PROJ, "raw", "geo")
REPORT = os.path.join(PROJ, "logs", "geo_pbmc_st6gal1_report.txt")
GENES = ["ST6GAL1", "CD19", "MS4A1", "CD79A", "CD22", "ACTB"]
ILMN = {"ILMN_1653120": "ST6GAL1", "ILMN_1756501": "ST6GAL1", "ILMN_2384496": "ST6GAL1"}
SETS = [("GSE73953", "GPL4133"), ("GSE125818", "GPL10558"), ("GSE58539", "GPL10558")]
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 200:
        return dest
    last = ""
    for i in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=180) as r, open(dest + ".part", "wb") as f:
                while True:
                    b = r.read(1 << 20)
                    if not b:
                        break
                    f.write(b)
            os.replace(dest + ".part", dest)
            return dest
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
            time.sleep(4 * (i + 1))
    log(f"    下载失败：{last[:200]}  {url}")
    return None


def stub(acc):
    return acc[:-3] + "nnn"


def platform_map(gpl):
    """返回 {探针ID: 基因名}（只留GENES里的）"""
    dest = os.path.join(GEO, f"{gpl}_table.txt")
    p = fetch(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=data", dest)
    mp = {}
    if not p:
        return mp
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        head = None
        for line in f:
            if line.startswith(("^", "!", "#")):
                continue
            v = line.rstrip("\n").split("\t")
            if head is None:
                head = [h.upper() for h in v]
                cand = [i for i, h in enumerate(head) if h in ("GENE_SYMBOL", "SYMBOL", "GENE SYMBOL", "ILMN_GENE")]
                if not cand:
                    log(f"    {gpl} 平台表找不到基因名列：{head[:15]}")
                    return mp
                gi = cand[0]
                continue
            if len(v) > gi and v[gi].strip().upper() in GENES:
                mp[v[0]] = v[gi].strip().upper()
    return mp


def welch(a, b):
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a) / (len(a) - 1)
    vb = sum((x - mb) ** 2 for x in b) / (len(b) - 1)
    se = math.sqrt(va / len(a) + vb / len(b))
    if se == 0:
        return float("nan")
    t = (ma - mb) / se
    return math.erfc(abs(t) / math.sqrt(2))      # 正态近似，样本小，仅供参考


def group(title, chars):
    s = (title + " " + chars).lower()
    if re.search(r"healthy|control|\bhc\b|normal|_c$|\bhs\b", s):
        return "对照"
    if re.search(r"membranous|\bmn\b", s):
        return "膜性肾病"
    if re.search(r"iga", s):
        return "IgAN"
    return "未知"


def do_array(acc, gpl):
    log("")
    log(f"--- {acc}（{gpl}）")
    p = fetch(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(acc)}/{acc}/matrix/{acc}_series_matrix.txt.gz",
              os.path.join(GEO, f"{acc}_series_matrix.txt.gz"))
    if not p:
        return
    mp = dict(ILMN) if gpl == "GPL10558" else {}
    mp.update(platform_map(gpl))
    log(f"    平台表里找到目标探针 {len(mp)} 个")
    titles, chars, ids, data, intable = [], [], [], {}, False
    with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            v = line.rstrip("\n").split("\t")
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in v[1:]]
            elif line.startswith(("!Sample_characteristics_ch1", "!Sample_source_name_ch1")):
                c = [x.strip('"') for x in v[1:]]
                chars = c if not chars else [a + "; " + b for a, b in zip(chars, c)]
            elif line.startswith(("!Series_summary", "!Sample_data_processing")) and len(_out) < 4000:
                if line.startswith("!Sample_data_processing"):
                    log("    数据处理：" + v[1].strip('"')[:200])
            elif line.startswith("!series_matrix_table_begin"):
                intable = True
            elif line.startswith("!series_matrix_table_end"):
                intable = False
            elif intable:
                if not ids:
                    ids = [x.strip('"') for x in v[1:]]
                    continue
                pid = v[0].strip('"')
                if pid in mp:
                    try:
                        data[pid] = [float(x) if x not in ("", "null", "NA") else float("nan") for x in v[1:]]
                    except ValueError:
                        pass
    if not chars:
        chars = [""] * len(titles)
    grp = [group(t, c) for t, c in zip(titles, chars)]
    cnt = {g: grp.count(g) for g in sorted(set(grp))}
    log(f"    样本 {len(titles)} 个，分组：{cnt}")
    for t, c, g in list(zip(titles, chars, grp))[:6]:
        log(f"      例：[{g}] {t[:60]} | {c[:80]}")
    if not data:
        log("    没有取到目标探针的数据")
        return
    for pid in sorted(data, key=lambda k: (mp[k], k)):
        vals = data[pid]
        by = {}
        for x, g in zip(vals, grp):
            if not math.isnan(x):
                by.setdefault(g, []).append(x)
        txt = "  ".join(f"{g} {sum(a) / len(a):.3f}（n={len(a)}）" for g, a in sorted(by.items()))
        p_ = welch(by.get("IgAN", []), by.get("对照", []))
        d = (sum(by["IgAN"]) / len(by["IgAN"]) - sum(by["对照"]) / len(by["对照"])) if by.get("IgAN") and by.get("对照") else float("nan")
        log(f"    {mp[pid]:<8} {pid:<14} {txt}   IgAN-对照 差={d:+.3f}  P≈{p_:.2g}")


def do_listing(acc):
    log("")
    log(f"--- {acc} 附件清单（单细胞，不下载）")
    for name in ("filelist.txt", ""):
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(acc)}/{acc}/suppl/{name}"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=60) as r:
                txt = r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            log(f"    {url} 取不到：{type(e).__name__}: {e}")
            continue
        if name:
            for ln in txt.splitlines()[:80]:
                log("    " + ln[:200])
        else:
            for m in re.finditer(r'href="([^"?/][^"]*)"[^\n]*?(\d{4}-\d\d-\d\d \d\d:\d\d)\s+([\d.]+[KMG]?)', txt):
                log(f"    {m.group(1)}  {m.group(3)}")
        return


def main():
    os.makedirs(GEO, exist_ok=True)
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    log("A. 芯片数据（数值为各数据集自带的标准化表达值，多为log2；差值为IgAN减对照）")
    for acc, gpl in SETS:
        try:
            do_array(acc, gpl)
        except Exception as e:  # noqa: BLE001
            log(f"!! {acc} 异常：{type(e).__name__}: {e}")
            for ln in traceback.format_exc().splitlines()[-5:]:
                log("    " + ln)
    log("")
    log("B. 单细胞数据")
    for acc in ("GSE285335", "GSE189481"):
        do_listing(acc)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
