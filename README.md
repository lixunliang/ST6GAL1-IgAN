# ST6GAL1 and IgA nephropathy: analysis code

Code for "Genetic dissection of the ST6GAL1 locus separates B-cell ST6GAL1 expression from IgG sialylation in susceptibility to IgA nephropathy".

All analyses use publicly available data: GWAS and QTL summary statistics, and de-identified single-cell RNA-sequencing count matrices from GEO (GSE285335) for one analysis. No individual-level genotype data are used.

## Setup

Scripts read from and write to a project directory given by the environment variable `MR_PROJ`; reference LD (1000 Genomes) and deCODE files are read from `MR_OLD_PROJ`. The defaults in the scripts are the authors' local paths and should be overridden, for example:

```bash
export MR_PROJ=/path/to/project          # Windows: set MR_PROJ=D:\project
export MR_OLD_PROJ=/path/to/reference
```

Software used: Python 3.11 with NumPy (1.26), pandas (2.x), SciPy (1.11), openpyxl (3.1) and Matplotlib (3.8); Node.js 20+ with the `docx` package (9.x) for building Word documents. Version numbers are those used by the authors; nearby versions are expected to work. Scripts that download data need internet access to the listed public servers.

## Run order

1. `00_build_olink_cis_instruments.py` (Olink instruments from Sun et al. 2023 ST9/ST16)
2. Downloads: `131`, `134`, `142`, `143`, `146`
3. Screen: `144` (SomaScan), `145` (Olink)
4. Colocalization of candidates: `147`
5. ST6GAL1 follow-up: `148` (serum IgA, reverse MR) → `149` (instruments, epitope checks, GTEx) → `150`–`152` (IgG N-glycome and two-signal analysis) → `153`–`154` (expression QTL) → `155` (PheWAS) → `156` (East Asian LD) → `157`–`159` (GEO and single-cell) → `160` (sensitivity analyses)
6. Figures and tables: `figures_tables/make_figs.py`, `make_figs2.py`, `make_supp.py`, then the `build_*.js` scripts

Scripts `98`, `116`, `122` and `141` contain shared functions imported by later scripts.

## Data to download (see Methods for sources)

| Directory under MR_PROJ | Content |
|---|---|
| raw/outcome/igan_kiryluk2023/ | IgA nephropathy GWAS (Kiryluk et al. 2023): combined, European, East Asian, without Immunochip |
| raw/outcome/serum_iga_liu2022/ | Serum IgA GWAS (Liu et al. 2022) |
| raw/outcome/igg_glycome_klaric2020/ | IgG N-glycome GWAS (Klarić et al. 2020; doi:10.7488/ds/2481) |
| raw/protein/ukbppp/ | UKB-PPP Olink summary statistics (tar files per protein) |
| raw/reference/ukbppp_sun2023/ | Sun et al. 2023 Supplementary Tables (41586_2023_6592_MOESM3_ESM.xlsx) |
| MR_OLD_PROJ/raw/reference/ld/ | 1000 Genomes phase 3 PLINK files (EUR, EAS) |
| MR_OLD_PROJ/extracted/cis/ | deCODE SomaScan cis regions |
| FinnGen R12 endpoints | downloaded by 143 |

## Scripts

| Script | Purpose |
|---|---|
| 00_build_olink_cis_instruments.py | Olink cis instruments from Sun et al. 2023 ST9/ST16 |
| 98_mr_step1_panukb.py | Shared MR functions (harmonisation, IVW, weighted median, MR-Egger) |
| 116, 122, 141 | Shared helpers (PheWeb/GTEx queries, Olink tar scanning) |
| 131, 134, 142, 143, 146 | Downloads and data probes |
| 144_proteome_to_kidney_diseases.py | SomaScan cis-MR and colocalization for IgA nephropathy |
| 145_olink_igan_mr.py | Olink cis-MR screen |
| 147_olink_igan_coloc.py | Regional colocalization for 24 candidate proteins |
| 148_st6gal1_chain.py | ST6GAL1 to serum IgA; reverse MR |
| 149_perk_st6gal1_instruments.py | LD with missense variants, instrument-level estimates, GTEx eQTL |
| 150–152 | IgG N-glycome: file probe, 77-trait analysis, two-signal conditional analysis |
| 153–154 | Expression QTL: GTEx and eQTL Catalogue (remote tabix) |
| 155 | Phenome-wide lookups |
| 156 | East Asian LD and conditional analysis |
| 157–159 | Public transcriptome data (GEO), single-cell pseudobulk |
| 160 | Colocalization prior sensitivity, conditional colocalization, rs7634389 |
| figures_tables/ | Figures, supplementary tables and Word documents |

Comments and console messages in the scripts are partly in Chinese.
