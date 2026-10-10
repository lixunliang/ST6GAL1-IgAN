# -*- coding: utf-8 -*-
"""155步：ST6GAL1两个信号的全表型扫描（FinnGen R10/R12、BBJ、UKB）
  S1 rs4686837（循环ST6GAL1 / B细胞表达 / IgA肾病）；S2 rs7621161（IgG唾液酸化）
看两个位点各自还关联哪些疾病和化验指标（免疫、肾脏、血液、甲状腺等）。需要联网，约5分钟。
"""
import importlib.util
import os
import re
import sys

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT = os.path.join(PROJ, "logs", "st6gal1_phewas_report.txt")
VARIANTS = [("rs4686837", "S1", ("3", 187021889), ("3", 186739677), ("G", "A")),
            ("rs7621161", "S2", ("3", 187009382), ("3", 186727170), ("C", "A"))]
KEY = re.compile(r"iga|nephr|glomer|kidney|renal|arthrit|lupus|scleros|thyroid|asthma|immunoglob|globulin|lymph|leuk|myeloma|"
                 r"protein|albumin|crohn|colitis|celiac|coeliac|psoria|diabetes|allerg|eosin|neutro|monocyt|platelet|"
                 r"white blood|sjogren|sarcoid|spondyl|hepat|liver|infect|tonsil", re.I)
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def pv(r):
    for k in ("pval", "pvalue", "pValue", "p"):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            continue
    return 1.0


def name(r):
    return str(r.get("phenostring") or r.get("phenotype") or r.get("phenocode") or r.get("pheno") or "?")[:75]


def main():
    spec = importlib.util.spec_from_file_location("_m116", os.path.join(HERE, "116_phewas_finemap_eqtl.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.log = lambda s="": None
    os.makedirs(m.PHE_DIR, exist_ok=True)
    log("beta方向对应各库的alt等位基因（见每行的 ref-alt）；原始JSON存到 raw\\phewas\\")
    for rs, tag, hg38, hg37, (a1, a2) in VARIANTS:
        log("")
        log(f"===== {tag} {rs}")
        for src, tmpl, build in m.PHEWEB:
            c, p = hg38 if build == 38 else hg37
            got = None
            for x, y in ((a1, a2), (a2, a1)):
                try:
                    obj = m.fetch_json(tmpl.format(c=c, p=p, a=x, b=y), os.path.join(m.PHE_DIR, f"{src}_{rs}_{x}_{y}.json"))
                except Exception as e:  # noqa: BLE001
                    obj = None
                    log(f"  {src}: {type(e).__name__}: {e}")
                recs = m.find_records(obj) if obj is not None else []
                if recs:
                    got = (x, y, recs)
                    break
            if not got:
                log(f"  {src:<12} 没取到结果")
                continue
            x, y, recs = got
            recs = sorted(recs, key=pv)
            thr = 0.05 / max(len(recs), 1)
            log(f"  {src:<12} ref-alt {x}-{y}，{len(recs)} 个表型（Bonferroni阈值 {thr:.1e}）；P最小的15个：")
            for r in recs[:15]:
                log(f"      P={pv(r):.2g}{'*' if pv(r) < thr else ' '} beta={r.get('beta', r.get('effect', ''))}  {name(r)}")
            sel = [r for r in recs[15:] if pv(r) < 0.01 and KEY.search(name(r))]
            if sel:
                log("    其余P<0.01且属免疫/肾脏/血液等关键词的：")
                for r in sel[:25]:
                    log(f"      P={pv(r):.2g}  beta={r.get('beta', r.get('effect', ''))}  {name(r)}")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
