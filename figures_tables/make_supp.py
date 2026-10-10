import csv, re, math, os
import pandas as pd, numpy as np
PROJ=os.environ.get("MR_PROJ", r"F:\孟德尔")
L=os.path.join(PROJ,"logs")+os.sep; R=os.path.join(PROJ,"results")+os.sep; O="supp/"
OUT={"IGAN_ALL":"IgAN combined","IGAN_EUR":"IgAN European","IGAN_ASN":"IgAN East Asian","IGAN_NOIC":"IgAN combined excl. Immunochip","IGAN_FG":"IgAN FinnGen R12"}
TR={"欧洲":"European","亚洲":"East Asian","合并":"Combined","去IC":"Combined excl. Immunochip","IgAN欧洲":"IgAN European","IgAN亚洲":"IgAN East Asian","IgAN合并":"IgAN combined","IgAN去IC":"IgAN combined excl. Immunochip","IgA水平":"Serum IgA",
    "B细胞":"B cells","记忆B":"Memory B cells","初始B":"Naive B cells","浆母细胞":"Plasmablasts","T细胞":"T cells","单核":"Monocytes","ST6GAL1蛋白":"Plasma ST6GAL1",
    "IgAN 欧洲":"IgAN European","IgAN 亚洲":"IgAN East Asian","IgAN 合并":"IgAN combined","IgAN 去IC":"IgAN combined excl. Immunochip"}
def tr(x):
    x=str(x)
    for k in sorted(TR,key=len,reverse=True): x=x.replace(k,TR[k])
    x=x.replace("细胞","cells").replace("肝细胞样","hepatocyte-like ").replace("肝","liver").replace("全血","whole blood").replace("脾","spleen").replace("甲状腺","thyroid").replace("小肠","small intestine").replace("单核","monocytes").replace("中性粒","neutrophils").replace("芯片","microarray").replace("浆母cells","plasmablasts").replace("组织","Bulk tissue").replace("髓系","Myeloid")
    return x
xl={}
# ---- SD1: Olink MR full
d=pd.read_csv(R+"olink_kidney_disease_mr.tsv",sep="\t"); d=d[d.outcome.isin(OUT)].copy(); d["outcome"]=d.outcome.map(OUT)
d=d.rename(columns={"pqtl_beta":"pQTL_beta","pqtl_log10p":"pQTL_-log10P","beta":"logOR_per_SD","bonferroni":"Bonferroni_sig","MHC":"in_MHC"}).sort_values(["outcome","p"])
xl["S1_Olink_MR"]=d
s=pd.read_csv(R+"kidney_disease_protein_mr.tsv",sep="\t"); s=s[s.outcome.isin(OUT)].copy(); s["outcome"]=s.outcome.map(OUT); xl["S2_SomaScan_MR"]=s.rename(columns={"beta":"logOR_per_SD"}).sort_values(["outcome","p"])
c=pd.read_csv(R+"kidney_disease_coloc.tsv",sep="\t"); c=c[c.outcome.isin(OUT)].copy(); c["outcome"]=c.outcome.map(OUT); xl["S2b_SomaScan_coloc"]=c
oc=pd.read_csv(R+"olink_igan_coloc.tsv",sep="\t"); oc["IgAN_set"]=oc.IgAN_set.map(tr); oc=oc.drop(columns=["note"]); xl["S1b_Olink_coloc"]=oc
pi=pd.read_csv(R+"perk_st6gal1_instruments.tsv",sep="\t"); pi["outcome"]=pi.outcome.map(tr); xl["S3_instruments"]=pi
eq=pd.read_csv(R+"st6gal1_eqtl_catalogue.tsv",sep="\t"); eq["dataset"]=eq.dataset.map(tr); eq["category"]=eq.category.map(tr); xl["S4_eQTL_Catalogue"]=eq
sc=pd.read_csv(R+"scrna_pbmc_st6gal1_pseudobulk.tsv",sep="\t"); sc["celltype"]=sc.celltype.map(tr); sc["group"]=sc.group.map({"H":"Healthy","E":"Early IgAN","L":"Late IgAN"}); xl["S6_scRNA_pseudobulk"]=sc
gl=pd.read_csv(R+"st6gal1_igg_glycome.tsv",sep="\t"); gl.columns=[tr(x).replace("PP4_glycan_IgAN_","PP4_glycan_x_") for x in gl.columns]; xl["S7_IgG_glycome"]=gl
ts=pd.read_csv(R+"st6gal1_two_signals.tsv",sep="\t"); ts["trait"]=ts.trait.map(tr); xl["S7b_two_signal_conditional"]=ts
ch=pd.read_csv(R+"st6gal1_chain.tsv",sep="\t"); ch["analysis"]=ch.analysis.map(lambda x: tr(x).replace("protein_to_Serum IgA(全部)","protein_to_serum_IgA_all").replace("protein_to_Serum IgA(不含deCODE)","protein_to_serum_IgA_excl_deCODE")); xl["S3b_IgA_reverse_upstream"]=ch
# PheWAS from report
rows=[]; sig=None; src=None
for ln in open(L+"st6gal1_phewas_report.txt",encoding="utf-8"):
    m=re.match(r"===== (S\d) (rs\d+)",ln)
    if m: sig=m.group(1)+" "+m.group(2); continue
    m=re.match(r"\s+(\S+)\s+ref-alt (\w)-(\w)",ln)
    if m: src=m.group(1); ea=m.group(3); continue
    m=re.match(r"\s+P=([\d.e+-]+)\*?\s+beta=(\S*)\s+(.*)",ln)
    if m and sig: rows.append([sig,src,ea,m.group(3).strip(),m.group(2),float(m.group(1))])
ph=pd.DataFrame(rows,columns=["signal","biobank","effect_allele","phenotype","beta","p"]).drop_duplicates(); xl["S5_PheWAS"]=ph
DESC={'S1_Olink_MR': 'Olink cis-pQTL Mendelian randomization for IgA nephropathy, all testable proteins (Wald ratio per SD protein)', 'S2_SomaScan_MR': 'SomaScan (deCODE) cis-pQTL Mendelian randomization for IgA nephropathy', 'S2b_SomaScan_coloc': 'SomaScan colocalization results', 'S1b_Olink_coloc': 'Olink regional colocalization for the 24 candidate proteins', 'S3_instruments': 'Instrument-level estimates for ST6GAL1 and EIF2AK3', 'S4_eQTL_Catalogue': 'eQTL Catalogue lookups of the two ST6GAL1 signals (beta for A allele)', 'S6_scRNA_pseudobulk': 'Pseudobulk expression per donor and cell type, GSE285335', 'S7_IgG_glycome': 'ST6GAL1 lead pQTL against 77 IgG N-glycan traits', 'S7b_two_signal_conditional': 'Marginal and conditional z for the two ST6GAL1 signals', 'S3b_IgA_reverse_upstream': 'Serum IgA, reverse MR and exposure analyses', 'S5_PheWAS': 'Phenome-wide lookups (top associations and immune/kidney-related hits)'}
with pd.ExcelWriter(O+"Supplementary_Data_1.xlsx") as w:
    pd.DataFrame({"sheet":list(xl),"content":[DESC[k] for k in xl]}).to_excel(w,sheet_name="README",index=False)
    for k,v in xl.items(): v.to_excel(w,sheet_name=k[:31],index=False)
print({k:len(v) for k,v in xl.items()})
# ---- compact tables for Word (CSV)
def f(x,n=2):
    return "" if pd.isna(x) else f"{x:.{n}f}"
def fp(p):
    if pd.isna(p): return ""
    if p<0.001:
        m,e=f"{p:.1e}".split("e"); return f"{m} × 10^{int(e)}"
    return f"{p:.3f}" if p<0.1 else f"{p:.2f}"
T={}
a=d[(d.outcome=="IgAN combined")&(d.FDR<0.05)].copy()
piv={}
for o in ["IgAN European","IgAN East Asian","IgAN FinnGen R12"]:
    piv[o]=d[d.outcome==o].set_index("Gene")[["logOR_per_SD","p"]]
pp4={(r.Gene,r.IgAN_set):r.PP4 for r in oc.itertuples()}
t1=[]
for r in a.itertuples():
    row=[r.Gene, r.rsid, "Yes" if r.in_MHC else "No", f(r.logOR_per_SD), fp(r.p), fp(r.FDR)]
    for o in ["IgAN European","IgAN East Asian","IgAN FinnGen R12"]:
        x=piv[o].loc[r.Gene] if r.Gene in piv[o].index else None
        if x is not None and isinstance(x,pd.DataFrame): x=x.iloc[0]
        row.append(f"{f(x.logOR_per_SD)} ({fp(x.p)})" if x is not None else "")
    row.append(f(pp4.get((r.Gene,"Combined")),2) if (r.Gene,"Combined") in pp4 else "")
    t1.append(row)
T["Supplementary Table 1. Plasma proteins associated with IgA nephropathy at FDR < 0.05 in the combined analysis (Olink cis-pQTL Wald ratio)"]=(["Protein","Lead pQTL","MHC","Log OR per SD (combined)","P","FDR","European log OR (P)","East Asian log OR (P)","FinnGen log OR (P)","PP4 (combined)"],t1,
 "Log OR, log odds ratio of IgA nephropathy per standard deviation higher plasma protein. PP4, posterior probability of a shared causal variant within 500 kb of the lead pQTL (computed for 25 prespecified candidates only). Full results for all 1,749 testable proteins are in Supplementary Data 1.")
cc=c[(c.PP4>0.8)].sort_values("PP4",ascending=False)
T["Supplementary Table 2. SomaScan cis-pQTL proteins colocalizing with IgA nephropathy (PP4 > 0.8)"]=(["Outcome","Aptamer","Protein","Log OR per SD","P","PP3","PP4","MHC"],
 [[r.outcome,r.SeqId,r.Gene,f(r.mr_beta),fp(r.mr_p),f(r.PP3),f(r.PP4),"Yes" if r.MHC else "No"] for r in cc.itertuples()],"Wald ratio at the strongest cis variant of each aptamer (deCODE). Full results in Supplementary Data 1.")
T["Supplementary Table 3. Instrument-level estimates for ST6GAL1 and EIF2AK3"]=(["Protein","Outcome","Instrument","Source","Protein beta","Wald estimate","SE","P"],
 [[r.Gene,r.outcome,r.instrument,r.source.replace("ST9_lead","Lead pQTL").replace("ST16_cond","Conditional pQTL").replace("pruned_r2_0.1","IVW, r² < 0.1"),f(float(r.protein_beta),3) if str(r.protein_beta)!="nan" else "",f(r.wald,3),f(r.se,3),fp(r.p)] for r in pi.itertuples()],
 "Wald estimates are per standard deviation of plasma protein (log OR for IgA nephropathy; standard deviation units for serum IgA). IVW, fixed-effect inverse-variance weighted estimate of instruments pruned to r² < 0.1.")
t4=[]
order=[]
for (ds,mt),g in eq.groupby(["dataset","molecular_trait_id"],sort=False):
    g=g.set_index("signal")
    if "S1" not in g.index or "S2" not in g.index: continue
    z1=g.loc["S1","beta_A"]/g.loc["S1","se"]; z2=g.loc["S2","beta_A"]/g.loc["S2","se"]; rr=0.442; dd=math.sqrt(1-rr*rr)
    t4.append([g.loc["S1","category"],ds,str(int(g.loc["S1","n"])),mt.replace("ENSG00000073849","gene"),
        f"{g.loc['S1','beta_A']:+.3f} ({fp(g.loc['S1','p'])})", f"{g.loc['S2','beta_A']:+.3f} ({fp(g.loc['S2','p'])})", f"{(z1-rr*z2)/dd:+.1f}", f"{(z2-rr*z1)/dd:+.1f}"])
T["Supplementary Table 4. Expression QTL effects of the two ST6GAL1 signals on ST6GAL1 in 32 eQTL Catalogue datasets"]=(["Category","Dataset","n","Trait","Signal 1 beta (P)","Signal 2 beta (P)","Signal 1 z | signal 2","Signal 2 z | signal 1"],t4,
 "Beta is for the A allele of rs4686837 (signal 1) and rs7621161 (signal 2). Conditional z-scores use r = 0.442 (1000 Genomes European). For microarray datasets each probe is shown separately; 'gene' denotes gene-level quantification.")
ph5=ph[(ph.p<1e-3)|ph.phenotype.str.contains("IgA",case=False)].sort_values(["signal","p"])
T["Supplementary Table 5. Phenome-wide associations of the two ST6GAL1 signals (P < 0.001, and IgA nephropathy)"]=(["Signal","Biobank","Effect allele","Phenotype","Beta","P"],
 [[r.signal,r.biobank,r.effect_allele,r.phenotype,r.beta,fp(r.p)] for r in ph5.itertuples()],"Phenome-wide significance thresholds (Bonferroni across phenotypes in each resource): FinnGen 2.0 × 10^-5, Biobank Japan 1.6 × 10^-4, UK Biobank 3.5 × 10^-5. Full lookups in Supplementary Data 1.")
t6=[]
for ctp in ["B cells","Memory B cells","Naive B cells","T cells","Monocytes"]:
    for g in ["ST6GAL1","CD22","C1GALT1","C1GALT1C1","ACTB"]:
        x=sc[(sc.celltype==ctp)&(sc.gene==g)&(sc.n_cells>=20)]
        med=lambda grp: x[x.group==grp].cpm
        def mw(a,b):
            from scipy.stats import mannwhitneyu
            return mannwhitneyu(a,b).pvalue if len(a)>2 and len(b)>2 else float("nan")
        allp=x[x.group!="Healthy"].cpm
        t6.append([ctp,g]+[f"{med(k).median():.1f} (n={len(med(k))})" for k in ["Healthy","Early IgAN","Late IgAN"]]+[fp(mw(allp,med("Healthy")))])
T["Supplementary Table 6. Pseudobulk expression in peripheral blood mononuclear cells of healthy donors and patients with IgA nephropathy (GSE285335)"]=(["Cell type","Gene","Healthy","Early IgAN","Late IgAN","P (all IgAN vs healthy)"],t6,
 "Values are median counts per million across donors; n, number of donors with at least 20 cells of that type. P from two-sided Mann–Whitney U test with donors as the unit (exact, scipy).")
g7=gl.sort_values("p")
T["Supplementary Table 7. Association of the ST6GAL1 lead pQTL locus with 77 IgG N-glycan traits"]=(["Trait","Wald (per SD ST6GAL1)","P","PP3 protein–glycan","PP4 protein–glycan","Min P in window","Top glycan variant","PP4 glycan–IgAN (Eur/EAS/combined)"],
 [[r["trait"],f(r["wald"],3),fp(r["p"]),f(r["PP3"]),f(r["PP4"]),fp(r["glycan_min_p"]),r["top_glycan"],"/".join(f(r[c],2) for c in [x for x in gl.columns if x.startswith("PP4_glycan_x_")][:3])] for _,r in g7.iterrows()],
 "Wald estimates at the lead pQTL are confounded by linkage with signal 2 for sialylation traits and should not be read as causal effects (see Fig. 3). Trait codes follow the UPLC derived-trait scheme of Klarić et al.")
import json; json.dump({k:v for k,v in T.items()},open(O+"tables.json","w"),ensure_ascii=False,default=str)
print({k[:30]:len(v[1]) for k,v in T.items()})
