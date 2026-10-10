import re, math, csv
from statistics import NormalDist
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
PROJ=os.environ.get("MR_PROJ", r"F:\孟德尔")
L=os.path.join(PROJ,"logs")+os.sep
plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"xtick.labelsize":7.5,"ytick.labelsize":7.5,
 "axes.linewidth":0.8,"xtick.direction":"out","ytick.direction":"out","axes.spines.top":False,"axes.spines.right":False,
 "lines.linewidth":0.8,"pdf.fonttype":42})
C={"b":"#0072B2","o":"#D55E00","g":"#009E73","k":"#333333","y":"#E69F00"}
def panel(ax,l): ax.text(-0.02,1.04,l,transform=ax.transAxes,fontsize=10,fontweight="bold",ha="right",va="bottom")
def save(fig,name):
    for ext,kw in (("png",dict(dpi=300)),("pdf",{})): fig.savefig(f"{name}.{ext}",bbox_inches="tight",**kw)
    plt.close(fig)
nd=NormalDist()
# ---------- Fig 1: candidates
txt=open(L+"olink_igan_coloc_report.txt",encoding="utf-8").read()
genes=["TNFSF13","TNFRSF13B","CFH","ST6GAL1","EIF2AK3","ITGAM"]
sets=[("欧洲","European",C["b"],"o"),("亚洲","East Asian",C["o"],"s"),("合并","Combined",C["k"],"D")]
rows={}
sec=txt.split("C. 共定位")[1].split("D. 汇总")[0]
cur=None
for ln in sec.splitlines():
    m=re.match(r"^  (\w+)",ln)
    if m and not ln.startswith("    "): cur=m.group(1); continue
    m=re.match(r"^    (\S+)\s+位点.*PP4=([\d.]+).*在IgA肾病 P=([\d.e+-]+)，Wald=([+-][\d.]+)",ln)
    if m and cur in genes:
        p=float(m.group(3)); w=float(m.group(4)); z=abs(nd.inv_cdf(max(p,1e-300)/2)); se=abs(w)/z if z>0 else float("nan")
        rows[(cur,m.group(1))]=(w,se,float(m.group(2)))
fig,(a,b)=plt.subplots(1,2,figsize=(6.9,2.9),gridspec_kw={"width_ratios":[2.2,1],"wspace":0.08})
for i,g in enumerate(genes):
    for j,(k,lab,col,mk) in enumerate(sets):
        if (g,k) not in rows: continue
        w,se,pp=rows[(g,k)]; y=i+(j-1)*0.24
        a.errorbar(w,y,xerr=1.96*se,fmt=mk,color=col,ms=3.5,elinewidth=0.8,capsize=0,label=lab if i==0 else None)
        b.scatter(j,i,s=10+110*pp,color=col,marker=mk,linewidths=0)
        b.text(j,i+0.42,f"{pp:.2f}",ha="center",va="center",fontsize=6.5)
a.axvline(0,color="grey",lw=0.6,ls="--"); a.set_yticks(range(len(genes))); a.set_yticklabels(["APRIL (TNFSF13)","TACI (TNFRSF13B)","CFH","ST6GAL1","PERK (EIF2AK3)","ITGAM"]); a.invert_yaxis()
a.set_xlabel("log OR for IgA nephropathy per SD higher plasma protein (95% CI)")
a.legend(frameon=False,ncol=3,loc="lower center",bbox_to_anchor=(0.5,1.0),handletextpad=0.2,columnspacing=1.0)
b.set_ylim(a.get_ylim()); b.set_yticks([]); b.spines["left"].set_visible(False); b.set_xticks(range(3)); b.set_xticklabels(["EUR","EAS","All"]); b.set_xlim(-0.6,2.6)
b.set_xlabel("Colocalization PP4")
panel(a,"a"); panel(b,"b"); save(fig,"Fig2_candidates")
# ---------- Fig 2: two signals
t2=open(L+"st6gal1_two_signals_report.txt",encoding="utf-8").read()
def grab(name):
    for ln in t2.splitlines():
        if ln.strip().startswith(name):
            m=re.search(r"S1: beta=([+-][\d.]+) z=\s*([+-][\d.]+).*?S2: beta=([+-][\d.]+) z=\s*([+-][\d.]+)(?:.*?S1 z=\s*([+-][\d.]+).*?S2 z=\s*([+-][\d.]+))?",ln)
            return [float(x) if x else None for x in m.groups()]
traits=[("ST6GAL1蛋白","Plasma ST6GAL1"),("IGP29","IgG sialylation\n(FA2G2S1 ratio)"),("IGP24","IgG sialylation\n(fucosylated)"),("IGP54","Neutral FA2BG2"),("IgAN 欧洲","IgA nephropathy\n(European)")]
fig,(a,b)=plt.subplots(1,2,figsize=(6.9,2.8),gridspec_kw={"width_ratios":[1.7,1],"wspace":0.5})
import numpy as np
x=np.arange(len(traits)); wd=0.36
v=[grab(k) for k,_ in traits]
a.barh(x-wd/2,[abs(r[4]) for r in v],wd,color=C["b"],label="Signal 1 (rs4686837), conditioned on signal 2")
a.barh(x+wd/2,[abs(r[5]) for r in v],wd,color=C["o"],hatch="//",edgecolor="white",linewidth=0,label="Signal 2 (rs7621161), conditioned on signal 1")
a.axvline(1.96,color="grey",lw=0.6,ls="--"); a.set_xscale("symlog",linthresh=5); a.set_xlim(0,45); a.set_xticks([0,2,5,10,20,40]); a.set_xticklabels(["0","2","5","10","20","40"])
a.set_yticks(x); a.set_yticklabels([l for _,l in traits],fontsize=7.5); a.invert_yaxis(); a.set_xlabel("|Conditional z|")
a.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.45,1.08),fontsize=7)
ig=[("IgAN 欧洲","European"),("IgAN 亚洲","East Asian"),("IgAN 合并","Combined")]
for i,(k,lab) in enumerate(ig):
    r=grab(k)
    for j,bz in enumerate(((r[0],r[1]),(r[2],r[3]))):
        beta,z=bz; se=abs(beta/z); y=i+(j-0.5)*0.3
        b.errorbar(beta,y,xerr=1.96*se,fmt="o" if j==0 else "s",color=C["b"] if j==0 else C["o"],ms=3.5,elinewidth=0.8)
b.axvline(0,color="grey",lw=0.6,ls="--"); b.set_yticks(range(3)); b.set_yticklabels([l for _,l in ig]); b.invert_yaxis()
b.set_xlabel("log OR for IgA nephropathy\nper A allele (95% CI)")
panel(a,"a"); panel(b,"b"); save(fig,"Fig3_two_signals")
# ---------- Fig 3: cell-type eQTL
t3=open(L+"st6gal1_eqtl_catalogue_report.txt",encoding="utf-8").read().split("互相校正后z")[1]
cat=None; pts=[]
names={"B":"B cells / plasmablasts","LCL":"LCL","肝":"Liver","组织":"Bulk tissues","髓系":"Myeloid","T/NK":"T / NK cells"}
for ln in t3.splitlines():
    m=re.match(r"\s+-- (\S+)",ln)
    if m: cat=m.group(1); continue
    m=re.search(r"([+-][\d.]+) / ([+-][\d.]+)(?:\s+\[(\w+)\])?\s*$",ln)
    if m and cat: pts.append((cat,float(m.group(1)),float(m.group(2))))
fig,ax=plt.subplots(figsize=(3.35,3.0))
cols=[C["b"],C["g"],C["y"],C["k"],C["o"],"#CC79A7"]; mks=["o","s","^","D","v","P"]
for (k,lab),col,mk in zip(names.items(),cols,mks):
    q=[p for p in pts if p[0]==k]
    ax.scatter([p[1] for p in q],[p[2] for p in q],s=16,color=col,marker=mk,label=f"{lab} (n={len(q)})",linewidths=0,zorder=3)
for f in (ax.axvline,ax.axhline):
    f(0,color="grey",lw=0.6); f(-1.96,color="grey",lw=0.5,ls="--"); f(1.96,color="grey",lw=0.5,ls="--")
ax.set_xlabel("Conditional z, signal 1 (rs4686837-A)"); ax.set_ylabel("Conditional z, signal 2 (rs7621161-A)")
ax.set_xlim(-6,3); ax.set_ylim(-5,3); ax.legend(frameon=False,loc="upper left",bbox_to_anchor=(1.0,1.0),fontsize=7,handletextpad=0.2)
save(fig,"_old_celltype")
print(rows); print(v); print(len(pts))
