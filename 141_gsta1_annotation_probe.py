#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
141_gsta1_annotation_probe.py  —  GSTA1的组织定位、可成药性与变异层面的公开注释（2026-10-06）

目的：给"T2D → GSTA1 → eGFR"补上遗传学以外的证据层：
  1 Human Protein Atlas：GSTA1的组织/单细胞特异性、是否分泌入血、血中浓度
  2 GTEx v8：GSTA1及同家族GSTA2在各组织的表达中位数；rs2608629在肾髓质的eQTL
  3 Open Targets：靶点可成药性、已知药物、小鼠敲除表型、遗传约束；
     领头位点 6_52760996_G_A 所在的可信集关联了哪些性状（GWAS与各类QTL）
这些接口的返回结构我没法事先核实，所以这一步把原始JSON全部存下来，并打印摘要；
哪个接口失败会如实写明，不影响其他。

需要联网，一两分钟。
输出  raw\\annotation\\*.json ; logs\\gsta1_annotation_report.txt
"""
import json
import os
import sys
import time
import traceback
import urllib.error
import urllib.parse
import urllib.request

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
OUT = os.path.join(PROJ, "raw", "annotation")
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "gsta1_annotation_report.txt")
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-data-download"
GENES = {"GSTA1": ("ENSG00000243955", "ENSG00000243955.5"), "GSTA2": ("ENSG00000244067", "ENSG00000244067.2")}
VARIANT38 = "6_52760996_G_A"
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def http(url, data=None, tries=3):
    last = ""
    for i in range(1, tries + 1):
        try:
            h = {"User-Agent": UA, "Accept": "application/json"}
            if data is not None:
                h["Content-Type"] = "application/json"
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read().decode("utf-8", "replace")[:600]
            except Exception:  # noqa: BLE001
                pass
            last = f"HTTP {e.code} {body}"
            if e.code in (400, 403, 404, 422):
                break
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
        time.sleep(3 * i)
    log(f"    失败（{last[:500]}）：{url[:150]}")
    return None


def save(name, raw):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), "wb") as f:
        f.write(raw)
    try:
        return json.loads(raw.decode("utf-8", "replace"))
    except ValueError:
        return None


def brief(obj, n=1500):
    s = json.dumps(obj, ensure_ascii=False)
    return s[:n] + ("…" if len(s) > n else "")


def hpa():
    log("1. Human Protein Atlas")
    for g, (ens, _v) in GENES.items():
        raw = http(f"https://www.proteinatlas.org/{ens}.json")
        if raw is None:
            continue
        obj = save(f"hpa_{g}.json", raw)
        if not isinstance(obj, dict):
            log(f"  {g}：返回的不是JSON对象，已存原文")
            continue
        log(f"  {g}：共 {len(obj)} 个字段")
        for k, v in obj.items():
            kl = k.lower()
            if any(x in kl for x in ("specific", "secretome", "blood concentration", "subcellular", "protein class",
                                     "tissue expression cluster", "single cell", "distribution")):
                log(f"    {k}: {brief(v, 400)}")


def gtex():
    log("")
    log("2. GTEx v8")
    for g, (_e, ver) in GENES.items():
        raw = http("https://gtexportal.org/api/v2/expression/medianGeneExpression?"
                   + urllib.parse.urlencode({"gencodeId": ver, "datasetId": "gtex_v8", "itemsPerPage": 250}))
        if raw is None:
            continue
        obj = save(f"gtex_median_{g}.json", raw)
        rows = (obj or {}).get("data") or []
        rows = sorted(rows, key=lambda r: -float(r.get("median", 0) or 0))
        log(f"  {g} 表达中位数（TPM）最高的10个组织，共 {len(rows)} 个组织：")
        for r in rows[:10]:
            log(f"    {r.get('tissueSiteDetailId', ''):<40} {float(r.get('median', 0)):.1f}")
        for r in rows:
            if "Kidney" in str(r.get("tissueSiteDetailId", "")):
                log(f"    [肾] {r.get('tissueSiteDetailId')} {float(r.get('median', 0)):.1f}")
    for tis in ("Kidney_Medulla", "Kidney_Cortex", "Liver"):
        raw = http("https://gtexportal.org/api/v2/association/dyneqtl?" + urllib.parse.urlencode(
            {"gencodeId": GENES["GSTA2"][1], "variantId": "chr6_52760996_G_A_b38", "tissueSiteDetailId": tis,
             "datasetId": "gtex_v8"}))
        if raw is None:
            continue
        obj = save(f"gtex_eqtl_GSTA2_{tis}.json", raw) or {}
        log(f"  rs2608629 × GSTA2 在 {tis}: NES={obj.get('nes')}  P={obj.get('pValue')}")


def gql(name, query, variables):
    raw = http("https://api.platform.opentargets.org/api/v4/graphql",
               data=json.dumps({"query": query, "variables": variables}).encode("utf-8"))
    if raw is None:
        return None
    obj = save(name, raw)
    if obj and obj.get("errors"):
        log(f"    接口报错：{brief(obj['errors'], 600)}")
    return (obj or {}).get("data")


def opentargets():
    log("")
    log("3. Open Targets")
    q1 = """query($id:String!){ target(ensemblId:$id){ approvedSymbol approvedName biotype
      tractability{label modality value}
      geneticConstraint{constraintType oe oeLower oeUpper score}
      mousePhenotypes{modelPhenotypeLabel modelPhenotypeClasses{label}}
      knownDrugs{count rows{drug{name} phase mechanismOfAction disease{name}}}
      subcellularLocations{location} } }"""
    d = gql("opentargets_target_GSTA1.json", q1, {"id": GENES["GSTA1"][0]})
    t = (d or {}).get("target")
    if t:
        log(f"  {t.get('approvedSymbol')}：{t.get('approvedName')}")
        tr = [x for x in (t.get("tractability") or []) if x.get("value")]
        log("  可成药性为真的条目：" + ("；".join(f"{x['modality']}:{x['label']}" for x in tr) or "无"))
        log("  遗传约束：" + brief(t.get("geneticConstraint"), 500))
        mp = t.get("mousePhenotypes") or []
        log(f"  小鼠表型 {len(mp)} 条：" + "；".join(sorted({x.get('modelPhenotypeLabel', '') for x in mp})[:40]))
        kd = t.get("knownDrugs") or {}
        log(f"  已知药物记录 {kd.get('count')} 条：" + "；".join(sorted({(r.get('drug') or {}).get('name', '') + '(期' + str(r.get('phase')) + ')' for r in (kd.get('rows') or [])})[:30]))
    q2 = """query($id:String!){ variant(variantId:$id){ id rsIds mostSevereConsequence{label}
      credibleSets(page:{size:200,index:0}){ count rows{ studyLocusId studyType pValueMantissa pValueExponent beta
        study{ id traitFromSource nSamples projectId biosample{biosampleName} } } } } }"""
    d = gql("opentargets_variant_rs2608629.json", q2, {"id": VARIANT38})
    v = (d or {}).get("variant")
    if v:
        cs = v.get("credibleSets") or {}
        rows = cs.get("rows") or []
        log(f"  位点 {v.get('id')} {v.get('rsIds')}  后果：{(v.get('mostSevereConsequence') or {}).get('label')}；"
            f"所在可信集 {cs.get('count')} 个")

        def pv(r):
            try:
                return float(r["pValueMantissa"]) * 10 ** float(r["pValueExponent"])
            except (KeyError, TypeError, ValueError):
                return 1.0
        for typ in sorted({r.get("studyType") for r in rows}, key=str):
            sub = sorted([r for r in rows if r.get("studyType") == typ], key=pv)
            log(f"    [{typ}] {len(sub)} 个：")
            for r in sub[:25]:
                s = r.get("study") or {}
                bs = (s.get("biosample") or {}).get("biosampleName") or ""
                log(f"      P={pv(r):.1e} beta={r.get('beta')}  {str(s.get('traitFromSource'))[:60]}  {bs}  [{s.get('projectId') or s.get('id')}]")
    else:
        log("  位点查询没有返回结果（可能是变异ID写法或接口结构不同，原始返回已存）")


def main():
    try:
        sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
    except Exception:  # noqa: BLE001
        pass
    os.makedirs(LOG_DIR, exist_ok=True)
    for fn in (hpa, gtex, opentargets):
        try:
            fn()
        except Exception as e:  # noqa: BLE001
            log(f"!! 异常：{type(e).__name__}: {e}")
            for ln in traceback.format_exc().splitlines():
                log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")


if __name__ == "__main__":
    main()
