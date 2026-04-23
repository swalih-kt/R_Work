# Stacked Horizontal Bar Plot — Variant Counts per Gene by ACMG Classification

An R script that generates a publication-quality **stacked horizontal bar plot** showing the number of variants per gene, grouped and color-coded by ACMG pathogenicity classification. Genes are ordered by total variant count for easy visual ranking.

---

## Requirements

- **R v4.0 or later**
- R packages:
  - `ggplot2`
  - `dplyr`
  - `readr`

---

## Input File

**`Swalih_MD_plot_input - novel_bar.tsv`** — Tab-separated file with one row per variant.

| Column | Description |
|--------|-------------|
| `Gene` | Gene symbol (e.g. `BRCA1`, `TP53`) |
| `VCF` | Variant identifier (e.g. HGVS notation or rsID) |
| `ACMG` | ACMG pathogenicity classification (`Pathogenic`, `Likely pathogenic`, `VUS`, etc.) |
| `Database_Status` | Whether the variant is `Known` or `Novel` in public databases |

> **Note**: Column names are reassigned inside the script — the input file columns must be in the order above (Gene, VCF, ACMG, Database_Status) even if the header names differ.

---

## What the Script Does

### Step 1: Load and Rename Columns
Reads the TSV input and standardizes column names to `Gene`, `VCF`, `ACMG`, and `Database_Status` for consistency in downstream processing.

### Step 2: Summarize Variant Counts
Groups the data by **Gene × ACMG classification** and counts the number of variants in each combination. This produces one row per gene-ACMG pair with a `count` column.

### Step 3: Generate the Plot
Builds a stacked horizontal bar plot where:
- **Y-axis**: Genes, ordered by total variant count (highest at top)
- **X-axis**: Number of variants (range 0–28, ticks every 2)
- **Fill color**: ACMG classification category
- **Labels**: Variant count printed in white inside each stack segment
- **Border**: Black panel border for a clean, publication-ready appearance

### Step 4: Save Output
Exports the figure as a high-resolution PNG.

---

## Output

| File | Description |
|------|-------------|
| `Figure3.png` | Stacked horizontal bar chart — variant counts per gene by ACMG class |
| **Resolution** | 300 DPI |
| **Dimensions** | 2500 × 3000 px |

---

## Plot Design Notes

- Genes are sorted by **total variant count** using `reorder(Gene, count)` — the gene with the most variants appears at the top.
- Each bar segment is outlined in black (`color = "black"`) to visually separate ACMG categories within a gene.
- Count labels are centered inside each segment using `position_stack(vjust = 0.5)` and displayed in white for contrast.
- The x-axis is fixed to a maximum of 28 with breaks every 2 units — update `scale_x_continuous()` if your dataset has more variants per gene.
- `theme_minimal` is used as the base with a custom black panel border added for print-quality output.

---

## Notes

- The script summarizes at the **gene level** — if multiple variants share the same Gene × ACMG combination, they are counted together into a single bar segment.
- `Database_Status` (Known/Novel) is included in the input but not used in this plot — it can be used for additional filtering or faceting if needed.
- Update `setwd()` and the input filename to match your local file paths before running.

---

## References

- [ggplot2 Documentation](https://ggplot2.tidyverse.org/)
- [ACMG Variant Classification Guidelines](https://www.acmg.net/ACMG/Medical-Genetics-Practice-Resources/Practice-Guidelines/ACMG_Content.aspx?ItemNumber=7447)
