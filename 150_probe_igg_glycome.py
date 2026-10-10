# -*- coding: utf-8 -*-
"""150步：查看IgG N-糖组GWAS（Klaric 2020）的文件结构，只读不改"""
import glob
import gzip
import io
import os
import sys
import tarfile

try:
    sys.stdout.reconfigure(encoding="gbk", errors="replace", line_buffering=True)
except Exception:  # noqa: BLE001
    pass
PROJ = os.environ.get("MR_PROJ", r"F:\孟德尔")
D = os.path.join(PROJ, "raw", "outcome", "igg_glycome_klaric2020")
REPORT = os.path.join(PROJ, "logs", "igg_glycome_probe_report.txt")
_out = []


def log(s=""):
    print(s, flush=True)
    _out.append(str(s))


def head(fobj, name, n=6):
    raw = gzip.GzipFile(fileobj=fobj) if name.endswith(".gz") else fobj
    st = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    for _ in range(n):
        line = st.readline()
        if not line:
            break
        log("      " + line.rstrip("\n")[:300])


def main():
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    if not os.path.isdir(D):
        log(f"找不到文件夹 {D}")
    else:
        for p in sorted(glob.glob(os.path.join(D, "*"))):
            log(f"{os.path.basename(p)}  {os.path.getsize(p):,} 字节")
        for p in sorted(glob.glob(os.path.join(D, "*.txt"))):
            log(f"\n--- {os.path.basename(p)}")
            with open(p, "rb") as f:
                head(f, p, 40 if "readme" in p.lower() else 6)
        for p in sorted(glob.glob(os.path.join(D, "*.tar*"))):
            log(f"\n--- {os.path.basename(p)} 的内容（读完整个压缩包约需几分钟）")
            try:
                with tarfile.open(p, "r:*") as tf:
                    k = 0
                    for m in tf:
                        k += 1
                        log(f"  {m.name}  {m.size:,}")
                        if m.isfile() and k <= 3:
                            head(tf.extractfile(m), m.name)
                    log(f"  共 {k} 个条目")
            except Exception as e:  # noqa: BLE001
                log(f"  !! 打不开：{type(e).__name__}: {e}")
    with open(REPORT, "w", encoding="utf-8") as f:
        f.write("\n".join(_out) + "\n")
    log(f"\n报告：{REPORT}")


if __name__ == "__main__":
    main()
