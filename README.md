# DEGify

Automated Bulk RNA-Seq Differential Expression Analysis Framework

## Overview

DEGify is an automated R framework for bulk RNA-seq downstream analysis. 
It streamlines differential expression analysis, visualization, and 
functional enrichment from count matrices using DESeq2.

DEGify provides:

- Differential expression analysis
- Volcano plots
- PCA visualization
- Heatmaps
- Gene Ontology enrichment
- Automated result exporting
- Reproducibility summaries

## Installation

```r
# Install required packages
install.packages(c(
  "devtools"
))

# Install DEGify from GitHub
devtools::install_github("MohamedElsisii/DEGify")


Even before GitHub release,
this structure is good.

---

# Example Workflow ⭐

VERY important.

```md id="rgctx2906"
## Example Workflow

```r
library(DEGify)

# Load count matrix
counts <- read.csv(
  "counts.csv",
  row.names = 1
)

# Load metadata
metadata <- read.csv(
  "metadata.csv",
  row.names = 1
)

# Run DEGify
results <- run_degify(
  counts = counts,
  metadata = metadata,
  padj_cutoff = 0.05,
  logfc_cutoff = 1,
  ontology = "BP"
)


---

# Outputs ⭐

```md id="rgctx2907"
## Outputs

DEGify automatically generates:

- Volcano plots
- PCA plots
- Heatmaps
- GO enrichment plots
- DEG result tables
- Analysis summary files
- Reproducibility session information

## Generated Files

DEGify exports results into the selected output directory:

- All_DEGs.csv
- Significant_DEGs.csv
- Volcano_Plot.png
- PCA_Plot.png
- Heatmap.png
- GO_Enrichment.csv
- GO_Barplot.png
- GO_Dotplot.png
- Analysis_Summary.txt
- sessionInfo.txt

## Citation

If you use DEGify in your research, please cite the package appropriately.