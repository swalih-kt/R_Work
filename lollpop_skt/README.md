# HLCS Mutation Lollipop Plot

An R-based visualization pipeline for plotting protein-level mutations in the **HLCS** gene using a lollipop plot. Mutations are annotated with ACMG pathogenicity classification, mutation type (color-coded by RANK), and overlaid on protein domain architecture. Built using the `trackViewer` Bioconductor package.

---

## Requirements

- **R v4.0 or later**
- R packages:
  - `trackViewer`
  - `grid`

---

## Input File

**`HLCS_lollipop_Final.tsv`** — Tab-separated file containing curated HLCS variants with the following columns:

| Column | Description |
|--------|-------------|
| `HGVS_NM` | HGVS nucleotide notation (transcript-level, e.g. `NM_001352514.2:c.2611_2614del`) |
| `frequency` | Number of times the variant has been independently curated and classified under ACMG criteria from publicly available databases (e.g. ClinVar, HGMD, literature) |
| `ACMG` | ACMG pathogenicity classification (`Pathogenic`, `Likely pathogenic`) |
| `AACHANGE` | Amino acid change notation (e.g. `K871Gfs*26`, `H816Y`) |
| `POS` | Amino acid position on the protein (used as x-axis coordinate) |
| `MutationEffect` | Functional consequence (`frameshift deletion`, `frameshift insertion`, `nonsynonymous`, `stopgain`) |
| `RANK` | Integer rank (1–5) encoding mutation type — used for color assignment in the plot |

### RANK to Mutation Type Mapping

| RANK | Mutation Type | Color |
|------|--------------|-------|
| 1 | Frameshift deletion | `#ff4000` (red-orange) |
| 2 | Frameshift insertion | `#990099` (purple) |
| 3 | Frameshift substitution | `#ff00d6` (magenta) |
| 4 | Nonsynonymous | `#118ab2` (blue) |
| 5 | Stopgain | `#66c2a5` (teal) |

> **Note**: `frequency` reflects how many times a variant has been independently curated and assigned an ACMG classification from publicly available databases such as ClinVar, HGMD, and published literature — it is **not** population allele frequency. A higher frequency indicates stronger evidence through repeated independent curation. The variant with the **highest frequency** is the only one labeled directly on the plot to avoid overcrowding.

---

## Protein Domain Annotations

Two HLCS protein domains are drawn on the lollipop baseline:

| Domain | Amino Acid Range | Color |
|--------|-----------------|-------|
| BPL_LPL Catalytic Domain | 601 – 798 | `#66c2a5` (teal) |
| BPL C-Terminal Domain | 816 – 862 | `#fc8d62` (orange) |

The full protein length is **873 amino acids** (`NM_001352514.2`).

---

## Output

A lollipop plot rendered in the R graphics device showing:

- **X-axis**: Amino acid position (1–873)
- **Y-axis**: Variant frequency (curation count)
- **Lollipop color**: Mutation type (by RANK)
- **Lollipop size**: Fixed (`cex = 0.6`)
- **Label**: Amino acid change shown only for the highest-frequency variant
- **Domain bar**: Protein domains drawn at the base of the plot
- **Legends**: Mutation type (left) and protein domain (right)
- **Title**: Gene name (`HLCS`) and transcript (`NM_001352514.2`)

---

## Notes

- The plot is drawn using `lolliplot()` from `trackViewer` with `yaxis = FALSE` — axis tick labels are added manually using `grid.text()` for full layout control.
- Lollipop head borders are set to `gray80` for visual separation between overlapping points.
- Legends are placed manually using `grid` viewport coordinates (`legendGrob` + `editGrob`) to allow precise positioning independent of plot margins.
- The label rotation is set to 45° (`label.parameter.rot = 45`) for readability at dense positions.

---

## References

- [trackViewer Bioconductor Package](https://bioconductor.org/packages/release/bioc/html/trackViewer.html)
- [HLCS Gene — OMIM](https://www.omim.org/entry/609018)
- [ClinVar HLCS Variants](https://www.ncbi.nlm.nih.gov/clinvar/?term=HLCS)
