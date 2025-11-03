# Quick Start Guide

## Installation

```r
devtools::install_github("RauLabUNC/genePackRat")
```

## Minimal Example

```r
library(genePackRat)

# Load QTL results
loci <- read.csv("significant_qtls.csv")
scans <- readRDS("qtl_scans.rds")
thresholds <- readRDS("qtl_thresholds.rds")

# Generate packet for first locus
generateLocusPacket(
  locus_cluster = loci[1, ],
  input_path = "data/",
  scan_data = scans,
  threshold_data = thresholds
)
```

## What You Need

### Required Data Files

Place your data in this structure:

```
data/
└── processed/joinLoci/
    ├── relational_tables/
    │   ├── genes_mouse.csv
    │   ├── orthology.csv
    │   ├── associations.csv
    │   ├── mouseGenePhenotypes.csv
    │   └── traitLoci.csv
    ├── geneTables/
    │   └── multTrait_cis-eQTL_nrvmExp.csv
    └── bulk_exp/
        └── rna_expression.csv
```

See README for required columns in each file.

### Optional: CC Founder Variants

```
data/processed/joinLoci/relational_tables/ccVariants/
├── gene_mutations.csv
└── snp_mutations.csv
```

## Output Structure

Each packet creates:

```
results/locus_packets/locus_chr1_100-110Mb/
├── gene_info_cluster_chr1_100-110Mb.xlsx
├── zoomPlots/
│   └── locus_zoom_*.pdf
├── README_summary.txt
└── founder_snp_table_*.csv
```

## Custom Paths

Override any file path:

```r
generateLocusPacket(
  locus_cluster = my_locus,
  input_path = "data/",
  rna_info_file = "custom/path/expression.csv",
  founder_mutations_file = NULL  # Skip variants
)
```

## Filtering Utility

Use `filterGenes()` for quick data filtering:

```r
# Keep only GWAS hits
hits <- filterGenes(
  inputTable = all_genes,
  referenceTable = gwas_results,
  by = "gene_id",
  joinType = "semi"
)

# Add expression filter
expressed <- filterGenes(
  inputTable = hits,
  filters = list(
    makeFilter("biotype", "==", "protein_coding"),
    makeFilter("TPM", ">", 100)
  )
)
```

## Getting Help

```r
?generateLocusPacket
?generateLocusZoomPlot
?generateGeneInfoExcel
?filterGenes
```

## Development

See `.claude_notes/` for detailed development notes and design decisions.
