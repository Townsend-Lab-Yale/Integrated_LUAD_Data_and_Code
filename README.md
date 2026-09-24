# Integrated LUAD Data and Code

Integrated lung adenocarcinoma (LUAD) dataset and R pipeline for: (1) cancer
effect size estimation for somatic variants, (2) variant × tumor-suppressor-gene
epistasis testing, and (3) comparison of human variant selection strength with
tumor burden in genetically engineered mouse models (GEMMs).

## 1. Data

### 1.1 Cohort origin and sequencing type

The integrated dataset merges LUAD cohorts spanning whole-genome (WGS),
whole-exome (WES), and targeted/panel sequencing, downloaded from cBioPortal
datahub or the NCI Genomic Data Commons, plus the Yale cohort from
institutional records.

| Folder | `Source` code | Sequencing | Study | Institution/project | *n* |
|---|---|---|---|---|---|
| `data/genie_9/` | `Genie` | Targeted (multiple panels) | AACR Project GENIE Consortium, 2017 [G1] | AACR Project GENIE v9.0.0 | 5,997 |
| `data/luad_fm-ad/` | `FM-AD` | Targeted (FoundationOne CDx) | Hartmaier et al., 2017 [G2] | Foundation Medicine Adult Cancer Dataset (GDC `FM-AD`) | 1,206 |
| `data/lung_msk_2017/` | `MSK2017` | Targeted (IMPACT 341, 410) | Jordan et al., 2017 [G3] | MSKCC | 372 |
| `data/nsclc_pd1_msk_2018/` | `MSK2018` | Targeted (IMPACT 341, 410, 468) | Rizvi et al., 2018 [G4] | MSKCC | 186 |
| `data/luad_tcga/` | `TCGA` | WES | TCGA Research Network, 2014 [G5] | The Cancer Genome Atlas | 515 |
| `data/luad_oncosg_2020/` | `OncoSG` | WES | Chen et al., 2020 [G6] | OncoSG | 300 |
| `data/lung_nci_2022/` | `NCI` | WGS | Zhang et al., 2021 [G7] | NCI Sherlock-Lung | 187 |
| `data/luad_broad/` | `Broad` | WES, WGS, or both | Imielinski et al., 2012 [G8] | Broad Institute | 180 |
| `data/yale_luad/` | `Yale` | WES | Kadara et al., 2017 [G9] | Yale University | 108 |
| `data/luad_cptac_2020/` | `CPTAC` | WES | Gillette et al., 2020 [G10] | CPTAC | 108 |
| `data/nsclc_tracerx_2017/` | `TracerX` | WES | Jamal-Hanjani et al., 2017 [G11] | TRACERx | 61 |
| `data/luad_mskcc_2015/` | `MSK2015` | WES | Rizvi et al., 2015 [G12] | MSKCC | 10 |
| `data/luad_tsp/` | `TSP` | Targeted | Ding et al., 2008 [G13] | Tumor Sequencing Project | — ¹ |

¹ TSP is merged into `integrated_data/` (`Source == "TSP"`) but excluded from
the 9,230-sample analytic cohort above (insufficient smoking-history
annotation)

**References:**
[G1] AACR Project GENIE Consortium. *Cancer Discov.* 2017;7(8):818-831.
[G2] Hartmaier RJ et al. *Cancer Res.* 2017;77(9):2464-2475.
[G3] Jordan EJ et al. *Cancer Discov.* 2017;7(6):596-609.
[G4] Rizvi H et al. *J Clin Oncol.* 2018;36(7):633-641.
[G5] Cancer Genome Atlas Research Network. *Nature.* 2014;511(7511):543-550.
[G6] Chen J et al. *Nat Genet.* 2020;52(2):177-186.
[G7] Zhang T et al. *Nat Genet.* 2021;53(9):1348-1359.
[G8] Imielinski M et al. *Cell.* 2012;150(6):1107-1120.
[G9] Kadara H et al. *Ann Oncol.* 2017;28(1):75-82.
[G10] Gillette MA et al. *Cell.* 2020;182(1):200-225.e35.
[G11] Jamal-Hanjani M et al. *N Engl J Med.* 2017;376(22):2109-2121.
[G12] Rizvi NA et al. *Science.* 2015;348(6230):124-128.
[G13] Ding L et al. *Nature.* 2008;455(7216):1069-1075.


### 1.2 Supporting reference files (`data/`)

| File(s) | Format | Contents |
|---|---|---|
| `bed_files/*.bed` | BED | Capture intervals for the targeted panels (TSP, MSK-IMPACT 341/410/468, FoundationOne); used by `cancereffectsizeR` for coverage-based mutation-rate denominators. |
| `gene_panels/foundation_one.txt`, `tsp.txt`, `msk341.txt`, `msk410.txt`, `msk468.txt` | gene list | Hugo Symbols in each named panel. |
| `gene_panels/all_panel_genes.txt`, `all_panel_samples.txt` | text list | Union of genes/samples covered by any panel. |
| `gene_panels/genie_panel_genes.txt`, `genie_panels_used.txt` | CSV | Panel gene content and per-sample panel assignment for GENIE. |
| `gene_panels/msk2017_panels_used.txt`, `msk2018_panels_used.txt` | text | Per-sample panel assignment, MSK 2017/2018 cohorts. |
| `genie_9/genomic_information.txt` | TSV | Captured intervals per GENIE panel (`SEQ_ASSAY_ID`). |
| `genie_9/data_clinical_sample.txt` | cBioPortal format | GENIE sample-level metadata. |
| `nsclc_tracerx_2017/case_lists/*.txt` | cBioPortal case-list | Sample IDs per TRACERx subset ("all", "sequenced"). |
| `hg38ToHg19.over.chain` | UCSC liftOver chain | GRCh38→hg19 conversion for `TCGA`/`FM-AD` (the only two sources on GRCh38). |
| `genes_list.txt` | gene list | Driver genes analyzed in the manuscript. |

### 1.3 Integrated / merged datasets (`integrated_data/`)

Built by harmonizing the sources into one sample table and one mutation
table; direct inputs to `code/01_CES_variants_calculation_and_plot.R`. All
files are comma-separated; the unlabeled first column is a row index from the
merge and can be ignored.

**`merged_luad_clinical.txt`** — one row per sample.

| Column | Meaning |
|---|---|
| `Sample ID` | Unique sample identifier, harmonized across cohorts. |
| `Smoker` | `True`/`False`/blank (ever-smoker/never-smoker/not reported). |
| `Stage` | Tumor stage as reported by the source cohort; **not** harmonized across studies. |
| `Progression Free Survival (months)` | PFS follow-up, numeric. |
| `Treatment` | Treatment info as reported (available for a subset of cohorts). |
| `Overall Survival (months)` | OS follow-up, numeric. |
| `Vital Status` | Coded per source cohort (commonly cBioPortal-style: `NED`, `AWD`, `DOD`, `DUK`; `-M` suffix = with metastasis). |
| `Patient ID` | Unique patient identifier (one patient may have multiple samples). |
| `Overall Survival (Months)` | Second OS column retained from the merge; equivalent to `Overall Survival (months)` where both are populated. |

**`merged_luad_maf.txt`** — one row per somatic variant.

| Column | Meaning |
|---|---|
| `Sample ID` | Links to `merged_luad_clinical.txt`. |
| `Chromosome` | No `chr` prefix. GRCh37/hg19 for all sources except `TCGA`/`FM-AD` (GRCh38, lifted over in `code/01_...R` via `data/hg38ToHg19.over.chain`). |
| `Start_Position` | 1-based coordinate (same build note as above). |
| `Mutation` | `chr:pos Ref>Alt`. |
| `Reference_Allele`, `Tumor_Seq_Allele2` | Reference / mutant allele. |
| `Source` | Cohort code (see Table 1.1). |
| `Variant_Classification` | MAF-standard consequence category. |
| `Panel` | Assay/panel ID (blank for exome/genome-wide sequencing). |

**`merged_final.txt`** — `merged_luad_maf.txt` joined with the clinical columns
from `merged_luad_clinical.txt` by `Sample ID`.

**`merged_luad_maf_TransvertToGene.csv`** — variant-to-gene annotation from
`cancereffectsizeR`'s MAF-loading step.

| Column | Meaning |
|---|---|
| `Unique_Patient_Identifier` | Sample ID. |
| `Chromosome`, `Start_Position`, `Reference_Allele`, `Tumor_Allele` | Post-liftover coordinates/alleles (hg19). |
| `variant_type` | e.g. `snv`. |
| `variant_id` | `chr:pos_Ref>Alt`. |
| `genes`, `top_gene` | Overlapping gene(s); primary gene assigned. |
| `top_consequence` | `Gene_AAchange` (e.g. `OR8B3_R292S`). |
| `prelift_chr`, `prelift_start` | Pre-liftover coordinates (`TCGA`/`FM-AD` only). |
| `liftover_strand_flip` | `TRUE` if liftover required a strand flip. |

**`sample_sequencingType_source_info.csv`** — one row per sample.

| Column | Meaning |
|---|---|
| `Unique_Patient_Identifier` | Sample ID. |
| `coverage` | `genome` / `exome` / `targeted`. |
| `covered_regions` | Covered-region reference set for mutation-rate correction. |
| `sig_analysis_grp` | Internal batch index for mutational-signature analysis. |
| `maf_source` | Fine-grained assay/panel ID. |

### 1.4 Mouse GEMM comparative data (`data/*_final.csv`)

Tumor-burden measurements from a GEMM system, digitized from a published figure
with [WebPlotDigitizer](https://automeris.io/WebPlotDigitizer/), used for the
human-vs-mouse "Comparative framework" analysis in `code/`.

Source: Blair LM, Juan JM, Sebastian L, et al. Oncogenic context shapes the
fitness landscape of tumor suppression. *Nat Commun.* 2023;14:6422



### 1.5 Precomputed R objects (`code/R_data/`)

Checkpoints so downstream steps don't require recomputing from raw MAF data:
`load_maf_cesa_WES_TGS_WGS.rds` (loaded/QC-filtered `CESAnalysis`, all
cohorts), `cesa_smoking.rds` / `cesa_nonsmoking.rds` (post effect-size
estimation, by smoking status), `epistasis_variant_TSG_output.Rdata` (output of
`02_epistasis_variant_TSG.R`), `gencode.v38lift37.basic.annotation.gtf.Rdata`
(pre-parsed GENCODE v38lift37/GRCh37 gene annotation).

## 2. Code

### Installation

Tested on macOS, R 4.3.0.

```r
install.packages(c('ggplot2', 'data.table', 'dplyr', 'rtracklayer', 'stringr', 'ggpubr', 'patchwork'))
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@v2.10.2", dependencies = TRUE, force = TRUE)
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release", dependencies = TRUE, force = TRUE)
```

### Workflow

Run in order from within `code/` (scripts use relative paths `../data/`,
`../integrated_data/`):

1. **`01_CES_variants_calculation_and_plot.R`** — loads/QC-filters MAF data for
   all cohorts, lifts over `TCGA`/`FM-AD` to hg19, estimates cancer effect
   sizes (smokers vs. never-smokers), compares to mouse GEMM tumor burden, and
   generates figures.
2. **`02_epistasis_variant_TSG.R`** — tests variant × TSG epistasis and
   generates the human-vs-mouse epistasis figures.
3. **`03_mouseStats_humanSCC_stat.R`** — mouse tumor-burden summary statistics
   by genotype, written to `code/SuppleStatistics/`.

## 3. License

- **Code** (`code/*.R`): [GNU GPL v3.0](LICENSE.txt).
- **Data**: third-party files (Section 1.1) retain their original source's
  license/terms (cBioPortal, NCI GDC, or AACR Project GENIE). 

