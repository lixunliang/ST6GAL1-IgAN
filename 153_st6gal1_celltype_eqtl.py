# -*- coding: utf-8 -*-
"""153步：ST6GAL1两个信号在不同组织/细胞里的eQTL
  S1 = chr3:187,021,889（GRCh38；rs4686837，循环ST6GAL1与IgA肾病的信号）
  S2 = chr3:187,009,382（GRCh38；IgG唾液酸化的信号）
A. GTEx v8（血、LCL、脾、肝、回肠、唾液腺等）
B. eQTL Catalogue（B细胞、浆细胞、LCL、单核、T细胞、肝等分选细胞数据集）
需要联网。接口返回的原始内容存到 raw\\annotation\\，接口格式和预想不同时报告里会留下原文片段。
"""
import importlib.util
import json
import os
import re
import sys
import traceback
import urllib.parse

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
HERE = os.path.dirname(os.path.abspath(__file__))
LOG_DIR = os.path.join(PROJ, "logs")
REPORT = os.path.join(LOG_DIR, "st6gal1_celltype_eqtl_report.txt")
GENE = "ST6GAL1"
SNPS = {"S1": ("3", 187021889, ("G", "A")), "S2": ("3", 187009382, ("C", "A"))}   # (另一等位基因, 参考A等位基因)
GTEX_TISSUES = ["Whole_Blood", "Cells_EBV-transformed_lymphocytes", "Spleen", "Liver", "Small_Intestine_Terminal_Ileum",
                "Minor_Salivary_Gland", "Colon_Transverse", "Lung", "Kidney_Cortex", "Thyroid", "Adipose_Subcutaneous"]
KEY = re.compile(r"b.?cell|lcl|lymphoblast|plasma|liver|hepat|blood|monocyte|cd4|cd8|t.?cell|neutrophil|nk|macrophage|dendritic|tonsil|ileum|colon", re.I)
EQC = "https://www.ebi.ac.uk/eqtl/api/v2"
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def js(raw):
    try:
        return json.loads(raw.decode("utf-8", "replace"))
    except Exception:  # noqa: BLE001
        return None


def gtex(m, vids):
    log("A. GTEx v8（NES对应variantId最后一个等位基因；P<0.05的标*）")
    raw = m.http("https://gtexportal.org/api/v2/reference/gene?" + urllib.parse.urlencode(
        {"geneId": GENE, "gencodeVersion": "v26", "genomeBuild": "GRCh38/hg38"}))
    o = js(raw) if raw else None
    try:
        gid = o["data"][0]["gencodeId"]
    except Exception:  # noqa: BLE001
        log("  查不到gencodeId：" + (m.brief(o, 300) if o else "无返回"))
        return None
    log(f"  {GENE} gencodeId = {gid}")
    for k, (c, pos, (oa, ea)) in SNPS.items():
        raw = m.http("https://gtexportal.org/api/v2/dataset/variant?" + urllib.parse.urlencode(
            {"chromosome": "chr" + c, "pos": pos, "datasetId": "gtex_v8"}))
        o = js(raw) if raw else None
        cand = []
        try:
            for d in o["data"]:
                log(f"  {k} GTEx里的变异：{d.get('variantId')}  {d.get('snpId')}  MAF={d.get('maf01')}")
                cand.append(d["variantId"])
        except Exception:  # noqa: BLE001
            log(f"  {k} 变异查询没有结果：" + (m.brief(o, 300) if o else "无返回"))
        if not cand:
            cand = [f"chr{c}_{pos}_{oa}_{ea}_b38", f"chr{c}_{pos}_{ea}_{oa}_b38"]
        vids[k] = cand
        for tis in GTEX_TISSUES:
            got = False
            for vid in cand:
                raw = m.http("https://gtexportal.org/api/v2/association/dyneqtl?" + urllib.parse.urlencode(
                    {"gencodeId": gid, "variantId": vid, "tissueSiteDetailId": tis, "datasetId": "gtex_v8"}), tries=1)
                o = js(raw) if raw else None
                if o and o.get("nes") is not None:
                    p = float(o["pValue"])
                    log(f"    {k} {tis:<36} NES={float(o['nes']):+.3f}  P={p:.2e} {'*' if p < 0.05 else ''}  （{vid}）")
                    got = True
                    break
            if not got:
                log(f"    {k} {tis:<36} 没取到")
    return gid.split(".")[0]


def eqc(m, ens, vids):
    log("")
    log("B. eQTL Catalogue（基因表达ge数据集；beta对应variant最后一个等位基因）")
    ds, start = [], 0
    while True:
        raw = m.http(f"{EQC}/datasets/?" + urllib.parse.urlencode({"quant_method": "ge", "size": 1000, "start": start}))
        o = js(raw) if raw else None
        if not isinstance(o, list) or not o:
            if start == 0:
                log("  数据集列表取不到：" + (m.brief(o, 500) if o is not None else "无返回"))
            break
        ds += o
        if len(o) < 1000:
            break
        start += 1000
    if not ds:
        return
    m.save("eqtl_catalogue_datasets_ge.json", json.dumps(ds).encode())
    log(f"  ge数据集 {len(ds)} 个；第一个的字段：{m.brief(ds[0], 400)}")
    sel = [d for d in ds if KEY.search(" ".join(str(d.get(x, "")) for x in ("tissue_label", "sample_group", "condition_label", "study_label")))]
    log(f"  按关键词选中 {len(sel)} 个")
    res, shown_err = [], 0
    for d in sel:
        did = d.get("dataset_id")
        lab = f"{d.get('study_label', '')}|{d.get('tissue_label', '')}|{d.get('condition_label', '')}"
        for k in SNPS:
            hit = None
            for vid in vids[k]:
                v = vid.replace("_b38", "")
                raw = m.http(f"{EQC}/datasets/{did}/associations?" + urllib.parse.urlencode(
                    {"gene_id": ens, "variant": v, "size": 100}), tries=1)
                o = js(raw) if raw else None
                if isinstance(o, list) and o:
                    rows = [r for r in o if str(r.get("gene_id", ens)) == ens and str(r.get("variant", v)) == v]
                    if rows:
                        hit = (rows[0], v)
                        break
                elif o is not None and shown_err < 3 and not isinstance(o, list):
                    shown_err += 1
                    log(f"    [接口返回] {did} {v}: {m.brief(o, 300)}")
            if hit:
                r, v = hit
                try:
                    res.append((k, lab, did, float(r.get("beta")), float(r.get("pvalue")), r.get("se"), v, r.get("median_tpm")))
                except (TypeError, ValueError):
                    log(f"    {did} {k} 字段解析不了：{m.brief(r, 300)}")
    if not res:
        log("  一个结果也没取到（可能接口格式变了，见上面的原文片段）")
        return
    for k in SNPS:
        log(f"  --- {k}（按P排序）")
        for _k, lab, did, b, p, se, v, tpm in sorted([r for r in res if r[0] == k], key=lambda r: r[4]):
            log(f"    {lab:<60} {did}  beta={b:+.3f}  P={p:.2e} {'*' if p < 0.05 else ''}  TPM中位数={tpm}  （{v}）")
    with open(os.path.join(PROJ, "results", "st6gal1_celltype_eqtl.tsv"), "w", encoding="utf-8") as f:
        f.write("signal\tdataset\tdataset_id\tbeta\tp\tse\tvariant\tmedian_tpm\n")
        for r in res:
            f.write("\t".join(str(x) for x in (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7])) + "\n")


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(os.path.join(PROJ, "results"), exist_ok=True)
    try:
        m = load("_m141", "141_gsta1_annotation_probe.py")
        m.log = log
        vids = {}
        ens = gtex(m, vids) or "ENSG00000073849"
        for k, (c, pos, (oa, ea)) in SNPS.items():
            vids.setdefault(k, [f"chr{c}_{pos}_{oa}_{ea}_b38", f"chr{c}_{pos}_{ea}_{oa}_b38"])
        eqc(m, ens, vids)
    except Exception as e:  # noqa: BLE001
        log(f"!! 异常：{type(e).__name__}: {e}")
        for ln in traceback.format_exc().splitlines():
            log("    " + ln)
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
