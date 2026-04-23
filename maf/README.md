# gnomAD Population Frequency & ACMG Variant Analysis Pipeline

An R-based pipeline for comparing allele frequencies of variants of interest across global gnomAD populations and the IndiGen cohort. The pipeline performs statistical testing (Fisher's exact test with Bonferroni correction) and produces a publication-ready bubble plot showing population-specific variant enrichment with ACMG pathogenicity classification.

---

## Requirements

- **R v4.0 or later**
- R packages:
  - `dplyr`
  - `ggplot2`
  - `ggforce`
  - `patchwork`

---

## Input Files

| File | Description |
|------|-------------|
| `all_variant.txt` | gnomAD variant table with allele counts and allele numbers across all populations (tab-separated) |
| `var_indigen.txt` | Variants of interest with IndiGen allele counts (`IndiGen_AC`, `IndiGen_AN`), ACMG classification, and gene annotation |

---

## Pipeline Steps

### Step 1: Data Loading and Column Selection

**Purpose**: Reads the full gnomAD variant table and retains only relevant columns — variant ID, allele counts (`_ac`), and allele numbers (`_an`) — while excluding homozygous (`_hom`) and hemizygous (`_hemi`) counts which are not needed for frequency comparisons.

| | |
|---|---|
| **Input** | `all_variant.txt` |
| **Output** | Filtered in-memory data frame (`gnomad1`) |

> **Important**: Column selection uses pattern matching on column names. Ensure gnomAD column naming follows the standard format. The `chr` prefix is appended to variant IDs to match hg38 coordinate style.

---

### Step 2: Column Renaming and Population Labelling

**Purpose**: Standardizes column names by replacing gnomAD's verbose naming convention with short population codes used throughout the pipeline.

| gnomAD label | Renamed to |
|---|---|
| `joint_` prefix | *(removed)* |
| `afr` | `AFR` |
| `amr` | `AMR` |
| `asj` | `ASJ` |
| `eas` | `EAS` |
| `nfe` | `NFE` |
| `fin` | `FIN` |
| `mid` | `MID` |
| `sas` | `SAS` |
| `ami` | `AMI` |
| `remaining` | `OTH` |

---

### Step 3: Merge with Variants of Interest and Compute Absent Alleles

**Purpose**: Filters the gnomAD dataset to only the variants listed in `var_indigen.txt` by merging on variant ID. Then computes `_ab` (absent allele) columns for each population — the number of alleles **not** carrying the variant — required for constructing 2×2 contingency tables in statistical testing.

| | |
|---|---|
| **Input** | `gnomad1`, `var_indigen.txt` |
| **Output** | Structured data frame (`df`) with AC, AN, and AN_ab columns for IndiGen and all gnomAD populations |

> **Important**: Absent alleles are calculated as `AN − AC` for each population. These values form the second row of the Fisher's exact test contingency matrix.

---

### Step 4: Fisher's Exact Test with Bonferroni Correction

**Purpose**: For each variant, tests whether the allele frequency in the **IndiGen cohort** is significantly different from each gnomAD population using Fisher's exact test. P-values are corrected for multiple comparisons using the Bonferroni method.

| | |
|---|---|
| **Input** | AC and AN_ab columns for IndiGen vs. each population |
| **Output** | One Bonferroni-corrected p-value column per population (e.g., `gnomad_AFR_bonf`, `gnomad_SAS_bonf`) |
| **Populations tested** | IndiGen, ALL, AFR, AMR, ASJ, EAS, FIN, MID, NFE, AMI, SAS, OTH |

> **Important**: The test is performed row-wise using a 2×2 matrix of [IndiGen AC, IndiGen AN_ab; Population AC, Population AN_ab]. Rows or columns summing to zero are skipped and assigned `NA` to avoid errors.

---

### Step 5: Minor Allele Frequency (MAF) Calculation

**Purpose**: Computes MAF for the IndiGen cohort and all gnomAD populations as `AC / AN`. Missing AC values are replaced with 0 and missing AN values are replaced with the column maximum to avoid division errors.

| | |
|---|---|
| **Input** | Cleaned `filtered_data` with AC/AN columns |
| **Output** | New `_MAF` columns for each population (e.g., `IndiGen_MAF`, `AFR_MAF`, `SAS_MAF`) |

> **Important**: Long VCF identifiers (>32 characters) are trimmed and suffixed with the number of trimmed nucleotides (e.g., `_5nts`) to keep plot labels readable.

---

### Step 6: Long-Format Data Preparation and ACMG Annotation

**Purpose**: Reshapes the wide MAF and Bonferroni p-value tables into long format for ggplot2 compatibility. Merges with ACMG pathogenicity classifications and gene names from `var_indigen.txt`.

| | |
|---|---|
| **Input** | MAF columns, Bonferroni p-value columns, `var_indigen.txt` |
| **Output** | `acmg` — long-format data frame with columns: VCF, Population, MAF, f_pvalue, ACMG, gene |

> **Important**: Population levels are set in a fixed order (IndiGen → ALL → AFR → AMR → ASJ → EAS → FIN → MID → NFE → AMI → SAS → OTH) to ensure consistent axis ordering in the plot.

---

### Step 7: Generate Population Enrichment Plot (Figure 7)

**Purpose**: Produces a two-panel publication-ready bubble plot showing population-specific variant frequencies with ACMG classification and statistical significance highlighted.

**Top panel** (`main_plot`): Bubble plot where each point represents a variant–population combination.
- **X-axis**: Variant ID
- **Y-axis**: Population
- **Bubble size**: MAF
- **Bubble colour**: ACMG pathogenicity class
- **Red circle overlay**: Marks variants with Bonferroni-corrected p-value < 0.05 (statistically enriched in IndiGen vs. that population)

**Bottom panel** (`bottom_axis_plot`): Colour-coded gene blocks beneath the x-axis showing which gene each variant belongs to.

| | |
|---|---|
| **Input** | `acmg` long-format data, gene block coordinates |
| **Output** | `Figure7.png` (4000 × 2800 px, 300 DPI) |

> **Important**: The two panels are combined using `patchwork` with a height ratio of 0.5:3 (gene bar : bubble plot). The gene block bar is rendered separately to allow independent x-axis scaling and tight layout control.

---

## Output Files

| File | Description |
|------|-------------|
| `Figure7.png` | Population-specific variant enrichment bubble plot with ACMG classification and significance overlay |

---

## Plot Legend Guide

| Visual Element | Meaning |
|---|---|
| Bubble colour | ACMG pathogenicity class (Pathogenic, Likely Pathogenic, VUS, etc.) |
| Bubble size | Minor Allele Frequency (MAF) in that population |
| Red open circle | Bonferroni-corrected p-value < 0.05 vs. IndiGen |
| Bottom colour bar | Gene to which the variant belongs |

---

## Notes

- The pipeline is designed for a **case cohort (IndiGen)** compared against gnomAD global and sub-population frequencies. All statistical comparisons are IndiGen vs. each gnomAD population, not cross-population.
- Variants absent from gnomAD (AC = 0, AN = 0) are retained in the analysis but will produce `NA` p-values and zero MAF for the missing populations.
- ACMG classifications must be pre-assigned in `var_indigen.txt` before running the pipeline.
- `patchwork` is loaded inside the plot block — move it to the library section at the top if reusing this script as a module.

---

## References

- [gnomAD Browser](https://gnomad.broadinstitute.org/)
- [IndiGen Programme](https://www.igib.res.in/content/indigen)
- [ACMG Variant Classification Guidelines](https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/Practice-Guidelines/ACMG_Content.aspx?ItemNumber=7447)
- [ggplot2 Documentation](https://ggplot2.tidyverse.org/)
- [patchwork Documentation](https://patchwork.data-imaginist.com/)
