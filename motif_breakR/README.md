# motifbreakR Analysis Pipeline

A pipeline for predicting the impact of regulatory variants (SNPs) on transcription factor (TF) binding motifs using **motifbreakR** in R. The pipeline loads variants from a BED file, scores them against HOCOMOCO and JASPAR 2024 motif databases, computes empirical p-values, and produces filtered output tables for downstream interpretation.

---

## Requirements

- **R v4.0 or later**
- Bioconductor packages:
  - `motifbreakR`
  - `MotifDb`
  - `BSgenome.Hsapiens.UCSC.hg38`
  - `SNPlocs.Hsapiens.dbSNP155.GRCh38`
  - `GenomicRanges`
  - `BiocParallel`
  - `BSgenome`
- CRAN packages:
  - `dplyr`

---

## Input Files

| File | Description |
|------|-------------|
| `clean_fixed.bed` | BED file of regulatory variants to analyze (chr, start, end, alleles) |

---

## Pipeline Steps

### Step 1: Load SNPs from BED File

**Purpose**: Reads variants from the input BED file and maps them to dbSNP155 rsIDs where available. Converts the variants into a `GRanges` object required by motifbreakR.

| | |
|---|---|
| **Input** | `clean_fixed.bed` |
| **Output** | `snps_output.tsv` — loaded SNP details with genomic coordinates and allele information |

> **Important**: `check.unnamed.for.rsid = TRUE` attempts to assign rsIDs to unnamed variants by matching coordinates against dbSNP155. Variants that cannot be matched will retain their BED-based identifiers.

---

### Step 2: Run motifbreakR (Without P-values)

**Purpose**: Scores each SNP against TF binding motifs from **HOCOMOCO v11** (core A, B, C) and **JASPAR 2024** databases. For each SNP, computes how much the alternate allele alters the motif match score compared to the reference allele.

| | |
|---|---|
| **Input** | Loaded SNPs (`snps`), MotifDb PWM subset |
| **Output** | `motifbreakR_results_no_pval.tsv` — all motif disruption scores without p-values |
| **Scoring method** | Information content (`ic`) |
| **Threshold** | `1e-4` (minimum motif score to retain a hit) |
| **Background** | Uniform nucleotide frequencies (A=C=G=T=0.25) |

> **Important**: `filterp = TRUE` pre-filters weak hits before full scoring, significantly reducing runtime on large variant sets. The motif databases used cover human TFs comprehensively — HOCOMOCO for curated human TFs and JASPAR 2024 for a broader cross-species set.

---

### Step 3: Calculate P-values

**Purpose**: Computes empirical p-values for each motif disruption score by comparing the observed score change against a null distribution. Adds statistical confidence to the effect size estimates from Step 2.

| | |
|---|---|
| **Input** | `results` object from Step 2 |
| **Output** | `motifbreakR_results_with_pval.tsv` — full results including p-values and effect classifications |
| **Parallelization** | `MulticoreParam` using all available CPU cores minus one |

> **Important**: This is the most computationally intensive step. Running on an HPC or multi-core workstation is strongly recommended for large variant sets.

---

### Step 4: Filter Strong Effect Hits

**Purpose**: Applies a stringent dual filter retaining only SNPs where both the **allele effect size** and the **p-value effect** are classified as `"strong"`. These represent the highest-confidence motif-disrupting variants.

| | |
|---|---|
| **Input** | `results_with_pval` from Step 3 |
| **Output** | `motifbreakR_final_hits_strong.tsv` — high-confidence motif disruption hits |
| **Filter** | `effect == "strong"` AND `pvalueEffect == "strong"` |

> **Important**: The dual filter (effect + p-value) is more conservative than filtering on either criterion alone. This reduces false positives and focuses interpretation on variants most likely to have a functional regulatory impact.

---

### Step 5: Per-SNP Summary

**Purpose**: Summarizes the strong-effect hits at the SNP level. For each variant, reports how many TFs are disrupted, which TFs are affected, and the mean and maximum allele effect sizes.

| | |
|---|---|
| **Input** | `motifbreakR_final_hits_strong.tsv` |
| **Output** | `SNP_summary.tsv` — one row per SNP with TF count, TF names, and effect statistics |
| **Sorted by** | Number of disrupted TFs (descending) |

> **Important**: SNPs disrupting many TFs simultaneously may indicate variants in highly conserved regulatory regions or motif clusters, and are strong candidates for functional follow-up.

---

### Step 6: Per-TF Summary

**Purpose**: Summarizes the strong-effect hits at the transcription factor level. For each TF, reports how many distinct SNPs disrupt its binding motif and lists the associated variant IDs.

| | |
|---|---|
| **Input** | `motifbreakR_final_hits_strong.tsv` |
| **Output** | `TF_summary_new.tsv` — one row per TF with SNP count and SNP IDs |
| **Sorted by** | Number of affected SNPs (descending) |

> **Important**: TFs disrupted by many SNPs in the dataset may represent recurrently affected regulatory hubs — particularly relevant in disease cohort studies where convergent TF disruption may point to shared regulatory pathways.

---

### Step 7: Filter by P-value Strength Only

**Purpose**: Alternative filtering that retains hits with a strong p-value regardless of effect size classification. Provides a broader hit set for comparison against the stricter dual filter applied in Step 4.

| | |
|---|---|
| **Input** | `results_with_pval` from Step 3 |
| **Output** | `pvalue_strong_hits.tsv` — hits filtered by p-value strength only |
| **Filter** | `pvalueEffect == "strong"` |

> **Important**: Comparing this output against Step 4 results can reveal hits with statistically significant p-values but moderate effect sizes — useful for secondary candidate prioritization.

---

### Step 8: Visualize a Single Variant *(Optional)*

**Purpose**: Subsets the results for a specific SNP of interest, recalculates its p-value individually, and generates a motif disruption plot. The plot displays the reference vs. alternate allele motif logo side by side, highlighting how the SNP alters TF binding affinity at that position.

| | |
|---|---|
| **Input** | `results` object, target rsID |
| **Output** | Motif logo plot (rendered in R graphics window or saved manually) |
| **Filter** | `effect = "strong"` — shows only strong-effect motifs for that SNP |

> **Important**: This step is intended for manual inspection of individual variants of interest — for example, a top candidate SNP from Step 5 or a variant with known functional relevance. Replace `"rs1479475149"` with any rsID present in your results.

---

## Output Files

| File | Description |
|------|-------------|
| `snps_output.tsv` | Loaded SNP coordinates and allele details |
| `motifbreakR_results_no_pval.tsv` | All motif disruption scores (no p-value) |
| `motifbreakR_results_with_pval.tsv` | Full results with p-values and effect classifications |
| `motifbreakR_final_hits_strong.tsv` | High-confidence hits (strong effect + strong p-value) |
| `SNP_summary.tsv` | Per-SNP summary of disrupted TFs and effect sizes |
| `TF_summary_new.tsv` | Per-TF summary of affected SNPs |
| `pvalue_strong_hits.tsv` | Hits filtered by p-value strength only |

---

## Notes

- **Motif databases used**: HOCOMOCO v11 (core A, B, C) and JASPAR 2024 — providing broad coverage of human TF binding preferences.
- **Scoring method**: Information content (`ic`) weights positions in the motif by their information content, giving more importance to highly conserved positions.
- **Background frequencies**: Uniform (0.25 each) — can be adjusted to genome-wide GC content for more accurate scoring in GC-rich regions.
- **Parallelization**: P-value calculation uses `MulticoreParam` with `detectCores() - 1` workers. On SLURM clusters, set `workers` to match your `--cpus-per-task` allocation.

---

## References

- [motifbreakR Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/motifbreakR.html)
- [HOCOMOCO v11](https://hocomoco11.autosome.org/)
- [JASPAR 2024](https://jaspar.elixir.no/)
- [MotifDb Bioconductor Page](https://bioconductor.org/packages/release/bioc/html/MotifDb.html)
