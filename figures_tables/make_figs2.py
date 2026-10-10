import re, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
exec(open("make_figs.py").read().split("nd=NormalDist()")[0])   # shared style, C, panel, save
# ---------------- Fig 1: design
fig,ax=plt.subplots(figsize=(6.9,6.0)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,txt,fc="#FFFFFF",ec="#333333",bold=False,fs=7.5):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.25,rounding_size=1.2",fc=fc,ec=ec,lw=0.8))
    ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=fs,fontweight="bold" if bold else "normal",linespacing=1.25)
def line(xs,ys): ax.plot(xs,ys,color="#333333",lw=0.8)
def arrow(x,y0,y1): ax.annotate("",xy=(x,y1),xytext=(x,y0),arrowprops=dict(arrowstyle="-|>",lw=0.8,color="#333333",mutation_scale=7))
B="#E8F1F8"; Og="#FBEDE3"; G="#E6F4EF"
box(2,80,30,17,"Plasma cis-pQTLs\nUKB-PPP (Olink):\n1,929 proteins\ndeCODE (SomaScan):\n1,829 aptamers",B,fs=6.8)
box(35,80,30,17,"IgA nephropathy GWAS\nKiryluk 2023: 10,146 cases,\n28,751 controls (European,\nEast Asian, combined)\nFinnGen R12: 822 cases",B,fs=6.8)
box(68,80,30,17,"Follow-up data\nSerum IgA (n = 41,263)\nIgG N-glycome (n = 8,090)\nGTEx v8; eQTL Catalogue\nPheWeb biobanks; GSE285335",B,fs=6.8)
line([17,17],[80,77.5]); line([50,50],[80,77.5]); line([17,50],[77.5,77.5]); arrow(33.5,77.5,75.5)
box(14,64,39,11,"Proteome-wide cis-MR (Wald ratio)\n1,749 testable proteins\n28 at FDR < 0.05 (8 in MHC)",Og,fs=7)
arrow(33.5,64,59.5)
box(14,46,39,13,"Colocalization (±500 kb), epitope checks\nAPRIL, TACI, CFH confirmed\nEIF2AK3: coding haplotype\nST6GAL1: OR 0.57 per SD, PP4 0.99",Og,fs=7)
arrow(33.5,46,43.5)
box(14,33,39,10,"ST6GAL1 locus dissection",G,bold=True)
line([83,83],[80,38]); line([53,83],[38,38])
xs=[1,21,41,61,81]; labels=["Serum IgA\nno effect\n(PP4 0.01)","Reverse MR\nno effect\n(P = 0.36)","IgG glycome\ntwo independent\nsignals","Sorted-cell eQTL\nsignal 1 in\nB cells / LCL","Biobank Japan\nsignal 1 with\nserum globulin"]
line([33.5,33.5],[33,27]); line([xs[0]+9,xs[-1]+9],[27,27])
for x,l in zip(xs,labels): arrow(x+9,27,22.5); box(x,8,18,14,l,fs=7)
ax.text(50,1.5,"Public single-cell PBMC data (GSE285335): no case–control difference in B-cell ST6GAL1",ha="center",fontsize=7,style="italic")
save(fig,"Fig1_design")
# ---------------- Fig 4: a = cell-type scatter (re-use), b = conditional z bars
t3=open(L+"st6gal1_eqtl_catalogue_report.txt",encoding="utf-8").read().split("互相校正后z")[1]
cat=None; pts=[]
names={"B":"B cells / plasmablasts","LCL":"LCL","肝":"Liver","组织":"Bulk tissues","髓系":"Myeloid","T/NK":"T / NK cells"}
for ln in t3.splitlines():
    m=re.match(r"\s+-- (\S+)",ln)
    if m: cat=m.group(1); continue
    m=re.search(r"([+-][\d.]+) / ([+-][\d.]+)(?:\s+\[(\w+)\])?\s*$",ln)
    if m and cat: pts.append((cat,float(m.group(1)),float(m.group(2))))
fig,(a,b)=plt.subplots(1,2,figsize=(6.9,3.0),gridspec_kw={"width_ratios":[1.15,1],"wspace":0.55})
cols=[C["b"],C["g"],C["y"],C["k"],C["o"],"#CC79A7"]; mks=["o","s","^","D","v","P"]
for (k,lab),col,mk in zip(names.items(),cols,mks):
    q=[p for p in pts if p[0]==k]
    a.scatter([p[1] for p in q],[p[2] for p in q],s=24,alpha=0.9,color=col,marker=mk,label=f"{lab} ({len(q)})",linewidths=0,zorder=3)
for f in (a.axvline,a.axhline):
    f(0,color="grey",lw=0.6); f(-1.96,color="grey",lw=0.5,ls="--"); f(1.96,color="grey",lw=0.5,ls="--")
a.set_xlabel("Signal 1 (rs4686837-A), conditional z"); a.set_ylabel("Signal 2 (rs7621161-A), conditional z")
a.set_xlim(-6,3); a.set_ylim(-5,3); a.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=6.5,handletextpad=0.1,columnspacing=0.6)
tr=[("IgA nephropathy\n(East Asian)",4.1,-0.9),("Non-albumin\nprotein",-5.5,-0.7),("Albumin/globulin\nratio",5.0,0.7),("Type 2 diabetes",-3.2,-1.4),("Monocyte count",-2.5,-1.3),("White cell count",-1.0,-2.8)]
y=np.arange(len(tr)); w=0.36
b.barh(y-w/2,[t[1] for t in tr],w,color=C["b"],label="Signal 1 | signal 2")
b.barh(y+w/2,[t[2] for t in tr],w,color=C["o"],hatch="//",edgecolor="white",linewidth=0,label="Signal 2 | signal 1")
for v in (-1.96,1.96): b.axvline(v,color="grey",lw=0.5,ls="--")
b.axvline(0,color="grey",lw=0.6); b.set_yticks(y); b.set_yticklabels([t[0] for t in tr],fontsize=7); b.invert_yaxis(); b.set_xlabel("Conditional z (A allele; East Asian LD)")
b.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.4,1.0),fontsize=6.5,ncol=1)
panel(a,"a"); panel(b,"b"); save(fig,"Fig4_celltype_phewas")
