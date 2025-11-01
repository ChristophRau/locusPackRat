# genePackRat: Gene Packet-based Rating for GWAS Gene Prioritization

genePackRat is an R package that provides a flexible, relational framework for prioritizing candidate genes from GWAS/QTL studies. It creates standardized "gene packets" - comprehensive evidence summaries that help researchers systematically identify the most promising candidate genes for experimental validation.

## Problem Statement

GWAS in model organisms often produces large loci (10-100+ Mb) containing dozens to hundreds of genes. Researchers need to integrate multiple data sources - expression, variants, phenotypes, orthologs - to identify the most likely causal genes. genePackRat provides a toolkit for flexible data integration and filtering based on relational principles.

## Key Features

- **Relational Filtering**: Filter any data based on relationships between tables
- **Flexible Column Mapping**: Works with your data structure, not prescriptive schemas
- **Modular Pipeline**: Chain operations together for complex workflows
- **Multi-format Output**: Excel workbooks, markdown reports, and publication-ready plots
- **General Purpose**: Works with any tabular data, not just genes

## Installation

### Prerequisites
genePackRat requires R >= 4.3.0 and depends on three CRAN packages (data.table, dplyr, openxlsx).

### Install from GitHub
```r
# Install with automatic dependency handling
install.packages("devtools")
devtools::install_github("RauLabUNC/genePackRat")
```

### Install from Local Directory

```r
# Option 1: Automatic installation (from genePackRat directory)
# This handles dependencies automatically:
source("install.R")

# Option 2: Manual installation
# First install dependencies:
install.packages(c("data.table", "dplyr", "openxlsx"))
# Then install genePackRat:
install.packages("/path/to/genepackrat", repos = NULL, type = "source")
```

### Optional Dependencies
For additional functionality (plotting, testing, etc.):

```r
# Optional but recommended
install.packages(c("testthat", "knitr", "rmarkdown", "RColorBrewer"))
```

## Quick Start

```r
library(genePackRat)

# Load your data - any data frame works
genes <- read.csv("my_genes.csv")
expression <- read.csv("expression_data.csv")
gwasHits <- read.csv("gwas_significant_genes.csv")

# Relational filtering - keep only GWAS hits
candidateGenes <- filterGenes(
  inputTable = genes,
  referenceTable = gwasHits,
  by = "gene_id",  # or by = c("symbol" = "gene_name") for different column names
  joinType = "semi"  # keeps only matching genes
)

# Add expression data and filter for high expression
expressedCandidates <- filterGenes(
  inputTable = candidateGenes,
  referenceTable = expression,
  by = "gene_id",
  joinType = "left",
  filters = list(
    makeFilter("TPM_heart", ">", 100),
    makeFilter("biotype", "==", "protein_coding")
  )
)

# Build comprehensive gene table with annotations
geneTable <- buildGeneTable(
  expressedCandidates,
  orthologyDf = mouseHumanOrthologs,
  expressionDf = tissueExpression,
  eqtlDf = eqtlResults,
  variantsDf = codingVariants
)

# Extract phenotypes matching your research area
cardiacPhenotypes <- extractPhenotypes(
  phenotypeData,
  useSets = c("cardiac", "metabolic")
)

# Create Excel workbook for curation
createGeneWorkbook(
  geneTable,
  outputFile = "locus_genes.xlsx"
)

# Generate summary report
report <- generateLocusReport(
  locusInfo,
  geneTable,
  outputFormat = "markdown"
)

# Create LocusZoom-style plot
plotLocus(
  chr = 8,
  start = 28000000,
  end = 32000000,
  genesDf = geneTable,
  outputFile = "locus_plot.png"
)
```

## Core Modules

### 1. filterGenes.R - Relational Filtering System
- `filterGenes()`: Unified filtering with relational joins and custom criteria
- `makeFilter()`: Create reusable filter specifications
- Supports semi/anti/inner/left joins between any tables
- Works with any column names via flexible `by` parameter

### 2. joinTables.R - Relational Data Integration
- `buildGeneTable()`: Combine gene annotations from multiple sources
- `aggregateGenesByLocus()`: Group genes by genomic regions
- Flexible column mapping for different data sources

### 3. phenotype.R - Phenotype Extraction
- `extractPhenotypes()`: Keyword-based phenotype matching
- `scorePhenotypeRelevance()`: Rank phenotypes by relevance
- `summarizePhenotypesByGene()`: Aggregate phenotype data

### 4. excel.R - Excel Workbook Generation
- `createGeneWorkbook()`: Multi-sheet Excel files for curation
- Automatic formatting and conditional highlighting
- Hyperlinks to external databases

### 5. summary.R - Report Generation
- `generateLocusReport()`: Markdown/HTML summary reports
- `identifyTopCandidates()`: Automated candidate prioritization
- `generateRecommendations()`: Follow-up suggestions

### 6. plotting.R - Visualization
- `plotLocus()`: LocusZoom-style genomic region plots
- `plotLocusDashboard()`: Multi-panel summary visualizations
- Support for plotgardener and base R graphics

## Design Philosophy

genePackRat follows tidyverse design principles:
- **Modular functions** that do one thing well
- **Consistent interfaces** with sensible defaults
- **Flexible inputs** supporting various identifier systems
- **Well-documented** with examples and vignettes
- **No hard dependencies** on specific databases or organisms

## Relational Filtering Approach

genePackRat treats gene prioritization as a relational data problem:

```r
# Example: Find cardiac GWAS hits with high expression
cardiacCandidates <- genes %>%
  filterGenes(referenceTable = gwasHits, by = "gene_id", joinType = "semi") %>%
  filterGenes(referenceTable = expression, by = "gene_id", joinType = "left",
              filters = list(makeFilter("TPM_heart", ">", 100)))
```

### Join Types Explained
- **`semi`**: Keep rows that match reference (like `%in%`)
- **`anti`**: Keep rows that DON'T match reference (like `!%in%`)
- **`inner`**: Merge matching rows from both tables
- **`left`**: Keep all input rows, add reference columns

### Flexible Column Mapping
```r
# Different column names? No problem!
filtered <- filterGenes(
  inputTable = mouseGenes,
  referenceTable = humanOrthologs,
  by = c("mouse_symbol" = "gene"),  # Map columns
  joinType = "inner"
)
```

## Use Cases

### Primary Use Case: QTL Gene Prioritization
Researchers with QTL intervals from mouse (CC, DO, BXD) studies who need to:
- Integrate multiple evidence types
- Create team-shareable curation documents
- Generate publication figures

### Secondary Use Cases
- Human GWAS with broad LD blocks
- Rat or other model organism QTL mapping
- Any scenario requiring systematic gene prioritization
- General relational filtering of biological data

## Data Model

genePackRat uses a flexible relational data model that adapts to your data structure:

```
Your Gene List ─┬─> Reference Table 1 ─> Filter ─> Candidates
                │         (join by X)
                │
                ├─> Reference Table 2 ─> Filter ─> Refined
                │         (join by Y)
                │
                └─> Reference Table N ─> Filter ─> Final List
                          (join by Z)
```


## Complete Example Workflow

Here's a real-world example showing the power of relational filtering:

```r
library(genePackRat)
library(dplyr)  # for pipe operator

# Load your data from various sources
allGenes <- read.csv("all_genes.csv")
gwasResults <- read.csv("gwas_significant.csv")
rnaseq <- read.csv("heart_rnaseq.csv")
eqtlData <- read.csv("heart_eqtl.csv")
mouseKnockouts <- read.csv("mgi_phenotypes.csv")

# Step 1: Start with GWAS hits
gwasGenes <- filterGenes(
  inputTable = allGenes,
  referenceTable = gwasResults,
  by = "gene_id",
  joinType = "semi"  # Keep only GWAS hits
)
# Started with 25,000 genes -> 487 GWAS hits

# Step 2: Add expression data and filter for expressed genes
expressedHits <- filterGenes(
  inputTable = gwasGenes,
  referenceTable = rnaseq,
  by = "gene_id",
  joinType = "left",  # Add expression columns
  filters = list(
    makeFilter("TPM", ">", 10),  # Expressed in heart
    makeFilter("padj", "<", 0.05)  # Differentially expressed
  )
)
# 487 genes -> 73 expressed and differential

# Step 3: Check for eQTL support
withEqtl <- filterGenes(
  inputTable = expressedHits,
  referenceTable = eqtlData,
  by = c("gene_id" = "gene"),  # Different column names
  joinType = "left",
  filters = list(
    makeFilter("eqtl_pvalue", "<", 1e-5)
  )
)
# 73 genes -> 28 with eQTL support

# Step 4: Exclude genes with conflicting phenotypes
finalCandidates <- filterGenes(
  inputTable = withEqtl,
  referenceTable = mouseKnockouts,
  by = c("gene_symbol" = "mgi_symbol"),
  joinType = "anti",  # EXCLUDE matches
  filters = list(
    function(df) !grepl("lethal", df$phenotype)
  )
)
# 28 genes -> 19 final candidates

# Create report
createGeneWorkbook(finalCandidates, "candidates.xlsx")
```

## Contributing

We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.

## Citation

If you use genePackRat in your research, please cite:

> Gural B, Kimball T, Luu A, Rau CD. genePackRat: A Flexible Framework for Prioritizing Candidate Genes from Broad GWAS Intervals and other Gene-Level Studies. *In preparation* (2025).

## License

MIT License - see LICENSE file for details

## Support

- Report issues: [GitHub Issues](https://github.com/RauLabUNC/genePackRat/issues)
- Documentation: [Package Website](https://raulabunc.github.io/genePackRat)
- Contact: bgural@unc.edu