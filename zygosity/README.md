# Genotype Count Heatmaps (Figure 6)

An R-based pipeline for visualizing zygosity distributions of variants and genes across individuals using heatmaps. Produces a three-panel publication-ready figure showing variant-level and gene-level genotype counts (heterozygous, homozygous, and unique heterozygous/homozygous) across a cohort.

---

## Requirements

- **R v4.0 or later**
- R packages:
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `ggpubr`
  - `grid`

---

## Input Files

Three tab-separated input files are required:

| File | Description |
|------|-------------|
| `f6_var.txt` | Variant-level genotype counts — one row per variant per zygosity category |
| `f6_gene_hethom.txt` | Gene-level heterozygous and homozygous counts per gene |
| `f6_gene_uhet_hom.txt` | Gene-level unique heterozygous and homozygous counts per gene |

### Expected Columns (all three files)

| Column | Description |
|--------|-------------|
| `variant` / `gene` | Variant identifier or gene name (long names are auto-trimmed to 32 characters) |
| `zygosity` | Zygosity category (e.g. `heterozygous`, `homozygous`) — used as x-axis |
| `Individuals_count` | Number of individuals in the cohort carrying that variant/gene in that zygosity state |

---

## Pipeline Steps

### Step 1: Load Input Data

**Purpose**: Reads all three input tables into R. Each file represents a different level of genotype summarization — variant-level, gene-level het/hom, and gene-level unique counts.

> **Important**: Long variant or gene names exceeding 32 characters are automatically trimmed and suffixed with the number of removed characters (e.g. `ENST00000123456789_longname_5nts`) to keep heatmap axis labels readable.

---

### Step 2: Generate Heatmaps

**Purpose**: A reusable `create_heatmap()` function builds a ggplot2 tile heatmap for each input file. Tile color intensity represents the number of individuals, and count labels are printed inside non-zero tiles.

Three heatmaps are generated:

| Panel | Input File | X-axis | Y-axis | Title |
|-------|-----------|--------|--------|-------|
| A | `f6_var.txt` | Zygosity | Variant | Variant-level Genotype Counts |
| B | `f6_gene_hethom.txt` | Zygosity | Gene | Gene-level Het/Hom Counts |
| C | `f6_gene_uhet_hom.txt` | Zygosity | Gene | Gene-level Unique Counts |

**Heatmap color scale**:
- Low count → `#f3eded` (near white)
- High count → `#4e9fe5` (blue)
- Missing / zero values → `#f3eded` (shown as blank tiles)

> **Important**: Tiles with a count of zero are left without a text label to reduce visual clutter — only non-zero cells display their count value.

---

### Step 3: Combine Panels

**Purpose**: Arranges the three heatmaps into a single composite figure using `ggarrange()`. Panel A occupies the left half; Panels B and C are placed side by side on the right half.

**Layout**:

```
┌──────────────┬──────────┬──────────┐
│              │          │          │
│      A       │    B     │    C     │
│  (Variant)   │ (Het/Hom)│ (Unique) │
│              │          │          │
└──────────────┴──────────┴──────────┘
```

Panel width ratio: A : (B+C) = 2 : 2

---

### Step 4: Save Output

**Purpose**: Exports the combined figure as a high-resolution PNG file suitable for publication.

| | |
|---|---|
| **Output file** | `papper.png` |
| **Resolution** | 300 DPI |
| **Dimensions** | 6000 × 4000 px |
| **Point size** | 18 |

---

## Output

| File | Description |
|------|-------------|
| `papper.png` | Three-panel genotype count heatmap (Figure 6) |

---

## Plot Design Notes

- **Tile borders** are drawn in grey (`color = "grey"`) with minimal size (`0.1`) to separate cells without visual dominance.
- **Axis labels** are rotated 45° on the x-axis for readability when zygosity categories are long.
- **Legend** is labeled `"Individuals count"` with an enlarged color bar (`1.2 cm height`) for clear interpretation in print.
- **Figure title** is left-aligned and bold, added via `annotate_figure()` after panel assembly.
- Long names are trimmed to a **maximum of 32 characters** using the `shorten_name()` helper function to prevent axis label overflow.

---

## Notes

- The script contains two versions of the pipeline — an early draft (red color scale) and the final version (blue color scale, larger font sizes). Only the **second/final version** should be used for publication figures.
- `aes_string()` is used for compatibility with older ggplot2 versions. If using ggplot2 v3.4+, this can be updated to use `.data[[col]]` syntax.
- Panel labels (`"A"`, `"B"`, `"C"`) are added automatically by `ggarrange()` with `font.label = list(size = 14)`.

---

## References

- [ggplot2 Documentation](https://ggplot2.tidyverse.org/)
- [ggpubr Documentation](https://rpkgs.datanovia.com/ggpubr/)
