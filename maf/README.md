
### Figure 7: Population-Specific Enrichment of Pathogenic Variants

This script processes gnomAD variant data and IndiGen cohort data to visualize population-specific enrichment of pathogenic variants.

**Inputs:**
- `all_variant.txt` (gnomAD variant data)
- `var_indigen.txt` (IndiGen variant data with ACMG annotations)

**Outputs:**
- `Figure7.png` (combined ACMG scatter plot with gene blocks)

**Workflow:**
1. Data cleaning and column renaming
2. Calculation of allele counts (AC/AN) and absent columns
3. Fisher’s exact test with Bonferroni correction
4. Minor allele frequency (MAF) calculation
5. Merge with ACMG annotations
6. Visualization using ggplot2 + patchwork
