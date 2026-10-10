# ST6GAL1 and IgA nephropathy: analysis code

Code for "Genetic dissection of the ST6GAL1 locus separates B-cell ST6GAL1 expression from IgG sialylation in susceptibility to IgA nephropathy".

All analyses use publicly available summary statistics. Scripts are numbered in the order they were run. Each script reads from and writes to a project directory given by the environment variable `MR_PROJ` (default `F:\孟德尔`); reference LD and deCODE files are read from `MR_OLD_PROJ`. Python 3.10+ with NumPy and pandas; SciPy and openpyxl for some steps; Matplotlib for figures; Node.js with the `docx` package for building the Word documents.

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
