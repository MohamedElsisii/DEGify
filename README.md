# DEGify

Automated RNA-seq downstream analysis framework for differential expression, visualization, GO enrichment, and HTML reporting in R.

---

## Overview

DEGify is an automated R framework for bulk RNA-seq downstream analysis using DESeq2. The package streamlines differential expression analysis, visualization, functional enrichment, and reproducible report generation from count matrices through an end-to-end automated workflow.

DEGify provides:

- Differential expression analysis
- Volcano plots
- PCA visualization
- Heatmaps
- Gene Ontology enrichment
- Automated HTML reports
- Result exporting
- Reproducibility summaries

---

## Installation

```r
install.packages("devtools")

devtools::install_github(
  "MohamedElsisii/DEGify"
)

# Load package
library(DEGify)

```
---

## Input File Format

DEGify requires:

### Count Matrix

A CSV file containing raw gene expression counts:

- Rows represent genes
- Columns represent samples
- First column contains gene IDs

Example:

| GeneID | Sample1 | Sample2 | Sample3 |
|---|---|---|---|
| ENSG000001 | 120 | 135 | 98 |
| ENSG000002 | 45 | 39 | 50 |

---

### Metadata File

A CSV file describing sample conditions:

- Rows represent samples
- Must contain a column named `condition`

Example:

| Sample | condition |
|---|---|
| Sample1 | Tumor |
| Sample2 | Tumor |
| Sample3 | Normal |

---

## Example Workflow

```r
library(DEGify)

# Load count matrix CSV
counts <- read.csv(
  system.file(
    "extdata",
    "example_counts.csv",
    package = "DEGify"
  ),
  row.names = 1
)

# Load metadata CSV
metadata <- read.csv(
  system.file(
    "extdata",
    "example_metadata.csv",
    package = "DEGify"
  ),
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
```

---

## Output Files

DEGify automatically generates:

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
- DEGify_Report.html

---

## Example Outputs

### Volcano Plot

<img width="3000" height="2400" alt="Volcano_Plot" src="https://github.com/user-attachments/assets/087cf96a-f18f-4b80-a61b-e017df53ab5b" />

### PCA Plot

<img width="2400" height="1800" alt="PCA_Plot" src="https://github.com/user-attachments/assets/ddcf0c3e-65be-4bbc-b831-4138b096e551" />

### Heatmap

<img width="2400" height="3000" alt="Heatmap" src="https://github.com/user-attachments/assets/afa8d6b6-4e52-409d-b86d-4a40321ffbc8" />

### Generated Report
DEGify automatically generates a publication-style HTML report summarizing differential expression analysis, visualization, and functional enrichment results.
<img width="940" height="386" alt="Screenshot 2026-05-10 145022" src="https://github.com/user-attachments/assets/2936cbb4-eac8-42c1-b2f0-1c2265277d3a" />
<img width="943" height="381" alt="Screenshot 2026-05-10 145036" src="https://github.com/user-attachments/assets/12e74864-cc7e-4ae6-a438-60adf9655d46" />


---

## Dependencies

DEGify uses:

- DESeq2
- clusterProfiler
- ggplot2
- enrichplot
- pheatmap
- org.Hs.eg.db

---

DEGify is under active development, and additional enrichment and visualization features will be added in future releases.

## License

MIT License

---

## Citation

If you use DEGify in your research, please cite the package appropriately.
