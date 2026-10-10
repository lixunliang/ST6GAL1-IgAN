# -*- coding: utf-8 -*-
"""154步：eQTL Catalogue（接口已停用，改为直接按位置读FTP上的索引文件）
ST6GAL1两个信号在分选细胞里的eQTL：
  S1 = chr3:187,021,889 G>A（rs4686837；循环ST6GAL1与IgA肾病）
  S2 = chr3:187,009,382 C>A（rs7621161；IgG唾液酸化）
对每个数据集：下载.tbi索引（缓存），只取这两个位置附近的一小段数据。beta对应A等位基因。
同时给出两位点互相校正后的z（r=+0.442，1000G欧洲，近似）。
需要联网（ftp.ebi.ac.uk 的https），约20-40分钟。可重复运行，已完成的数据集会跳过。
"""
import gzip
import json
import math
import os
import struct
import sys
import time
import traceback
import urllib.request
import zlib

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass

PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
LOG_DIR = os.path.join(PROJ, "logs")
RES = os.path.join(PROJ, "results")
TBI_DIR = os.path.join(PROJ, "raw", "annotation", "eqtl_catalogue_tbi")
DONE = os.path.join(PROJ, "raw", "annotation", "eqtl_catalogue_st6gal1_rows.json")
REPORT = os.path.join(LOG_DIR, "st6gal1_eqtl_catalogue_report.txt")
BASE = os.environ.get("EQC_BASE", "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats")
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-data-download"
GENE = "ENSG00000073849"
CHROM, S1, S2, R = "3", 187021889, 187009382, 0.442
DATASETS = [
    # (study, dataset, 标签, 样本量, 类别)
    ("QTS000026", "QTD000474", "Schmiedel_2018 B细胞(naive)", 91, "B"),
    ("QTS000036", "QTD000585", "Randolph_2021 B细胞", 89, "B"),
    ("QTS000037", "QTD000597", "Perez_2022 B细胞", 163, "B"),
    ("QTS000038", "QTD000608", "OneK1K B_naive", 876, "B"),
    ("QTS000038", "QTD000606", "OneK1K B_intermediate", 618, "B"),
    ("QTS000038", "QTD000607", "OneK1K B_memory", 726, "B"),
    ("QTS000038", "QTD000623", "OneK1K 浆母细胞", 795, "B"),
    ("QTS000007", "QTD000073", "CEDAR B细胞CD19(芯片)", 262, "B"),
    ("QTS000009", "QTD000080", "Fairfax_2012 B细胞CD19(芯片)", 281, "B"),
    ("QTS000013", "QTD000110", "GEUVADIS LCL", 445, "LCL"),
    ("QTS000029", "QTD000539", "TwinsUK LCL", 418, "LCL"),
    ("QTS000012", "QTD000095", "GENCORD LCL", 190, "LCL"),
    ("QTS000015", "QTD000221", "GTEx LCL", 147, "LCL"),
    ("QTS000015", "QTD000266", "GTEx 肝", 208, "肝"),
    ("QTS000023", "QTD000404", "PhLiPS 肝细胞样细胞", 79, "肝"),
    ("QTS000015", "QTD000326", "GTEx 脾", 227, "组织"),
    ("QTS000015", "QTD000356", "GTEx 全血", 670, "组织"),
    ("QTS000019", "QTD000373", "Lepik_2017 全血", 471, "组织"),
    ("QTS000015", "QTD000341", "GTEx 甲状腺", 574, "组织"),
    ("QTS000015", "QTD000321", "GTEx 小肠", 174, "组织"),
    ("QTS000038", "QTD000609", "OneK1K CD14单核", 560, "髓系"),
    ("QTS000002", "QTD000021", "BLUEPRINT 单核", 191, "髓系"),
    ("QTS000010", "QTD000081", "Fairfax_2014 单核(芯片)", 420, "髓系"),
    ("QTS000002", "QTD000026", "BLUEPRINT 中性粒", 196, "髓系"),
    ("QTS000038", "QTD000626", "OneK1K cDC2", 869, "髓系"),
    ("QTS000038", "QTD000612", "OneK1K CD4_Naive", 949, "T/NK"),
    ("QTS000038", "QTD000613", "OneK1K CD4_TCM", 954, "T/NK"),
    ("QTS000038", "QTD000615", "OneK1K CD8_Naive", 729, "T/NK"),
    ("QTS000038", "QTD000617", "OneK1K CD8_TEM", 946, "T/NK"),
    ("QTS000038", "QTD000620", "OneK1K NK", 951, "T/NK"),
    ("QTS000002", "QTD000031", "BLUEPRINT CD4 T", 167, "T/NK"),
    ("QTS000026", "QTD000479", "Schmiedel_2018 CD4 T(naive)", 88, "T/NK"),
]
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def pz(z):
    return math.erfc(abs(z) / math.sqrt(2))


def fetch(url, rng=None, tries=3):
    last = ""
    for i in range(1, tries + 1):
        try:
            if url.startswith("file:"):
                with open(url[5:], "rb") as f:
                    if rng:
                        f.seek(rng[0])
                        return f.read(rng[1] - rng[0] + 1)
                    return f.read()
            h = {"User-Agent": UA}
            if rng:
                h["Range"] = f"bytes={rng[0]}-{rng[1]}"
            with urllib.request.urlopen(urllib.request.Request(url, headers=h), timeout=120) as r:
                data = r.read()
                if rng and r.status != 206:
                    raise IOError(f"服务器不支持分段读取（状态 {r.status}）")
                return data
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
            time.sleep(3 * i)
    raise IOError(last)


def tbi_offset(tbi_raw, chrom, pos):
    """返回覆盖pos的16 kb窗口第一条记录的虚拟偏移量"""
    d = gzip.decompress(tbi_raw)
    if d[:4] != b"TBI\x01":
        raise ValueError("不是tabix索引")
    n_ref, _fmt, _cs, _cb, _ce, _meta, _skip, l_nm = struct.unpack_from("<8i", d, 4)
    names = d[36:36 + l_nm].split(b"\x00")[:n_ref]
    p = 36 + l_nm
    want = [chrom.encode(), b"chr" + chrom.encode()]
    for i in range(n_ref):
        n_bin, = struct.unpack_from("<i", d, p); p += 4
        for _ in range(n_bin):
            _b, n_chunk = struct.unpack_from("<Ii", d, p); p += 8 + 16 * n_chunk
        n_intv, = struct.unpack_from("<i", d, p); p += 4
        if names[i] in want:
            w = (pos - 1) >> 14
            offs = struct.unpack_from(f"<{n_intv}Q", d, p)
            for k in range(w, n_intv):
                if offs[k]:
                    return offs[k]
            return None
        p += 8 * n_intv
    return None


def bgzf_blocks(buf):
    """从buf里依次解出完整的BGZF块，返回 [(解压后的字节, 块长度)]"""
    out, p = [], 0
    while p + 18 <= len(buf):
        if buf[p:p + 4] != b"\x1f\x8b\x08\x04":
            raise ValueError("BGZF块头不对")
        xlen, = struct.unpack_from("<H", buf, p + 10)
        q, bs = p + 12, None
        while q < p + 12 + xlen:
            si1, si2, slen = buf[q], buf[q + 1], struct.unpack_from("<H", buf, q + 2)[0]
            if si1 == 66 and si2 == 67:
                bs = struct.unpack_from("<H", buf, q + 4)[0] + 1
            q += 4 + slen
        if bs is None:
            raise ValueError("BGZF块没有长度字段")
        if p + bs > len(buf):
            break
        out.append((zlib.decompress(buf[p + 12 + xlen:p + bs - 8], -15), bs))
        p += bs
    return out


def query(url, lo, hi):
    head = bgzf_blocks(fetch(url, (0, 70000)))[0][0].split(b"\n", 1)[0].decode().split("\t")
    ix = {h: i for i, h in enumerate(head)}
    tbi_path = os.path.join(TBI_DIR, os.path.basename(url) + ".tbi")
    if os.path.exists(tbi_path) and os.path.getsize(tbi_path) > 1000:
        with open(tbi_path, "rb") as f:
            tbi = f.read()
    else:
        tbi = fetch(url + ".tbi")
        os.makedirs(TBI_DIR, exist_ok=True)
        with open(tbi_path, "wb") as f:
            f.write(tbi)
    voff = tbi_offset(tbi, CHROM, lo)
    if voff is None:
        return head, [], len(tbi), 0
    coff, uoff = voff >> 16, voff & 0xFFFF
    rows, tail, first, nbytes, done = [], b"", True, 0, False
    step = 1 << 19
    while not done and nbytes < (1 << 27):
        buf = fetch(url, (coff, coff + step - 1))
        if not buf:
            break
        blocks = bgzf_blocks(buf)
        if not blocks:
            step *= 2
            continue
        for data, bs in blocks:
            coff += bs
            nbytes += bs
            if first:
                data, first = data[uoff:], False
            lines = (tail + data).split(b"\n")
            tail = lines.pop()
            for ln in lines:
                v = ln.decode("utf-8", "replace").split("\t")
                if len(v) < len(head):
                    continue
                if v[ix["chromosome"]].replace("chr", "") != CHROM:
                    done = True
                    break
                pos = int(v[ix["position"]])
                if pos > hi:
                    done = True
                    break
                if pos in (S1, S2) and v[ix["gene_id"]] == GENE:
                    rows.append({h: v[i] for h, i in ix.items()})
            if done:
                break
        if len(buf) < step:
            break
        time.sleep(1)
    return head, rows, len(tbi), nbytes


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    os.makedirs(os.path.dirname(DONE), exist_ok=True)
    cache = {}
    if os.path.exists(DONE):
        try:
            with open(DONE, "r", encoding="utf-8") as f:
                cache = json.load(f)
        except ValueError:
            cache = {}
    lim = int(os.environ.get("MR_LIMIT", "0") or 0)
    fails = 0
    for n, (qts, qtd, lab, ss, cat) in enumerate(DATASETS, 1):
        if lim and n > lim:
            break
        if qtd in cache:
            continue
        url = f"{BASE}/{qts}/{qtd}/{qtd}.all.tsv.gz"
        t0 = time.time()
        try:
            head, rows, tb, nb = query(url, min(S1, S2), max(S1, S2))
            cache[qtd] = rows
            log(f"[{n}/{len(DATASETS)}] {lab:<32} 取到 {len(rows)} 行（索引 {tb / 1e6:.1f} MB，数据 {nb / 1e6:.2f} MB，{time.time() - t0:.0f}s）")
            with open(DONE, "w", encoding="utf-8") as f:
                json.dump(cache, f)
            fails = 0
        except Exception as e:  # noqa: BLE001
            fails += 1
            log(f"[{n}/{len(DATASETS)}] {lab:<32} 失败：{type(e).__name__}: {str(e)[:200]}")
            if n == 1:
                for ln in traceback.format_exc().splitlines()[-6:]:
                    log("      " + ln)
            if fails >= 4:
                log("连续4个数据集失败，停止（多半是网络或服务器限制），可稍后重跑")
                break
        time.sleep(2)

    log("")
    log("结果（beta对应A等位基因；A在S1上降低循环ST6GAL1并升高IgA肾病风险，在S2上降低IgG唾液酸化）")
    log(f"  {'数据集':<32} n     S1 beta (P)             S2 beta (P)             互相校正后z：S1 / S2")
    out = []
    last = None
    for qts, qtd, lab, ss, cat in DATASETS:
        rows = cache.get(qtd)
        if rows is None:
            continue
        if cat != last:
            log(f"  -- {cat}")
            last = cat
        traits = sorted(set(r.get("molecular_trait_id", "") for r in rows))
        if not traits:
            log(f"  {lab:<32} {ss:<5} 没有ST6GAL1在这两个位点的记录（基因未表达或未检测）")
        for t in traits:
            z, txt = {}, ""
            for k, pos in (("S1", S1), ("S2", S2)):
                rr = [r for r in rows if r.get("molecular_trait_id") == t and int(r["position"]) == pos and r.get("alt") == "A"]
                if not rr:
                    txt += f"  {k} 没有              "
                    continue
                try:
                    b, se, p = float(rr[0]["beta"]), float(rr[0]["se"]), float(rr[0]["pvalue"])
                except ValueError:
                    txt += f"  {k} 无效              "
                    continue
                z[k] = b / se if se > 0 else 0.0
                txt += f"  {b:+.3f} ({p:.1e}){'*' if p < 0.05 else ' '}      "
                out.append([lab, cat, ss, t, k, b, se, p, rr[0].get("median_tpm", ""), rr[0].get("maf", "")])
            c = ""
            if len(z) == 2:
                d = math.sqrt(1 - R * R)
                c1, c2 = (z["S1"] - R * z["S2"]) / d, (z["S2"] - R * z["S1"]) / d
                c = f"{c1:+.1f} / {c2:+.1f}"
            suffix = f"  [{t}]" if len(traits) > 1 else ""
            log(f"  {lab:<32} {ss:<5}{txt}{c}{suffix}")
    with open(os.path.join(RES, "st6gal1_eqtl_catalogue.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tcategory\tn\tmolecular_trait_id\tsignal\tbeta_A\tse\tp\tmedian_tpm\tmaf\n")
        for r in out:
            f.write("\t".join(str(x) for x in r) + "\n")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
