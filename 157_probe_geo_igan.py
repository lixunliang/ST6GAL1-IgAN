# -*- coding: utf-8 -*-
"""157步：在GEO里查IgA肾病的人类表达数据集（找外周血/B细胞/扁桃体/单细胞的），只查目录不下载数据。联网，约2分钟。"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
REPORT = os.path.join(PROJ, "logs", "geo_igan_probe_report.txt")
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TERM = '("IgA nephropathy"[All Fields] OR "IgAN"[All Fields] OR "immunoglobulin A nephropathy"[All Fields] OR "IgA vasculitis"[All Fields]) AND "Homo sapiens"[Organism] AND gse[Entry Type]'
KEY = re.compile(r"blood|pbmc|b.?cell|lymphocyte|leukocyte|monocyte|tonsil|single.?cell|scrna|plasma cell|mononuclear|bone marrow|intestin|ileum|gut", re.I)
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def get(url):
    last = ""
    for i in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
            time.sleep(3 * (i + 1))
    log(f"  失败：{last[:200]}  {url[:120]}")
    return None


def main():
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    o = get(EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "gds", "term": TERM, "retmax": 500, "retmode": "json"}))
    ids = (o or {}).get("esearchresult", {}).get("idlist", [])
    log(f"GEO里人类IgA肾病相关的GSE共 {len(ids)} 个")
    rows = []
    for i in range(0, len(ids), 100):
        s = get(EU + "esummary.fcgi?" + urllib.parse.urlencode({"db": "gds", "id": ",".join(ids[i:i + 100]), "retmode": "json"}))
        for uid, r in ((s or {}).get("result", {}) or {}).items():
            if uid == "uids" or not isinstance(r, dict):
                continue
            rows.append(r)
        time.sleep(1)
    rows.sort(key=lambda r: -int(r.get("n_samples") or 0))

    def show(r, full):
        log(f"  {r.get('accession')}  样本数 {r.get('n_samples')}  {r.get('gdstype')}  平台 GPL{r.get('gpl')}  发布 {r.get('pdat')}  PMID {','.join(str(x) for x in (r.get('pubmedids') or []))}")
        log(f"      {str(r.get('title'))[:160]}")
        if full:
            log(f"      {str(r.get('summary'))[:420]}")
            sm = [str(x.get('title')) for x in (r.get('samples') or [])[:6]]
            if sm:
                log(f"      样本举例：{' | '.join(sm)[:300]}")
    hit = [r for r in rows if KEY.search(str(r.get("title")) + " " + str(r.get("summary")))]
    log("")
    log(f"A. 标题或摘要提到 血/PBMC/B细胞/扁桃体/单细胞/肠 的 {len(hit)} 个（详细）")
    for r in hit:
        show(r, True)
    log("")
    log("B. 其余（多为肾组织），只列标题")
    for r in rows:
        if r not in hit:
            show(r, False)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
