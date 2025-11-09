# locusPackRat

A flexible R package for managing genomic analysis projects with persistent data storage and integrated visualization.

## Overview

locusPackRat provides a project-based workflow for genomic analyses. Instead of passing data between functions, you create a persistent project that stores your genes or genomic regions along with any supplementary data. This allows for:

- **Persistent storage**: Load data once, use many times
- **Flexible integration**: Add any type of supplementary data
- **Reproducible outputs**: Regenerate visualizations and reports with different parameters
- **Multi-species support**: Works with human and mouse genomes (hg38, hg19, mm39, mm10)

## Installation

```r
# Install from GitHub
devtools::install_github("RauLabUNC/locusPackRat")

# For visualization features, you'll also need plotgardener from Bioconductor:
BiocManager::install("plotgardener")
```

### Dependencies
- R >= 4.3.0
- data.table
- jsonlite
- openxlsx (for Excel output)
- plotgardener (optional, for LocusZoom plots)

## Quick Start

### 1. Initialize a Project

```r
library(locusPackRat)

# With gene list
genes <- data.frame(
  gene_symbol = c("BRCA1", "TP53", "EGFR", "MYC"),
  expression = c(100, 250, 50, 300)
)

initPackRat(
  data = genes,
  mode = "gene",
  species = "human",
  genome = "hg38",
  project_dir = "my_project"
)

# Or with genomic regions (e.g., QTL results)
regions <- data.frame(
  chr = c(5, 10, 12),
  start = c(10000000, 50000000, 75000000),
  end = c(15000000, 55000000, 80000000),
  trait = c("height", "weight", "BMI")
)

initPackRat(
  data = regions,
  mode = "region",
  species = "mouse",
  genome = "mm39",
  project_dir = "qtl_project"
)
```

### 2. Add Supplementary Data

```r
# Add expression data
expression_data <- read.csv("expression.csv")
addRatTable(
  data = expression_data,
  table_name = "expression",
  link_type = "gene",
  link_by = "gene_symbol",  # Auto-detected if NULL
  project_dir = "my_project"
)

# Add QTL scan results
scan_results <- read.csv("qtl_scans.csv")
addRatTable(
  data = scan_results,
  table_name = "qtl_scans",
  link_type = "region",
  project_dir = "qtl_project"
)

# View available tables
listPackRatTables("my_project")
```

### 3. Generate Outputs

```r
# Create Excel report with all data
generateGeneSheet(
  format = "excel",
  include_supplementary = TRUE,  # Include all tables
  output_file = "results/gene_report.xlsx",
  project_dir = "my_project"
)

# Create filtered CSV
generateGeneSheet(
  format = "csv",
  filter_expr = "expression > 100",
  include_supplementary = c("expression", "phenotypes"),  # Specific tables
  output_file = "results/high_expression.csv",
  project_dir = "my_project"
)

# Generate LocusZoom plot (requires plotgardener)
generateLocusZoomPlot_v2(
  region_id = "region_1",
  scan_table = "qtl_scans",
  project_dir = "qtl_project"
)
```

## Key Functions

### Project Management
- `initPackRat()` - Initialize a new project with genes or regions
- `addRatTable()` - Add supplementary data tables
- `listPackRatTables()` - List available supplementary tables

### Output Generation
- `generateGeneSheet()` - Create formatted Excel or CSV outputs
- `generateLocusZoomPlot_v2()` - Create LocusZoom-style visualization

## Project Structure

Projects are stored in a `.locusPackRat` directory:

```
my_project/
└── .locusPackRat/
    ├── input/
    │   ├── genes.csv         # Core gene/region data
    │   └── orthology.csv      # Cross-species orthologs
    ├── supplementary/
    │   ├── expression.csv     # Added via addRatTable()
    │   ├── phenotypes.csv
    │   └── qtl_scans.csv
    ├── output/
    │   ├── gene_sheet.xlsx
    │   └── plots/
    └── config.json            # Project metadata
```

## Advanced Features

### Excel with Multiple Tabs

```r
generateGeneSheet(
  format = "excel",
  split_by = "criteria",
  split_criteria = list(
    "High_Expression" = "expression > 100",
    "Disease_Associated" = "!is.na(disease)",
    "Significant" = "p_value < 0.05"
  ),
  highlight_genes = c("BRCA1", "TP53"),
  project_dir = "my_project"
)
```

### Batch Processing

```r
# Generate plots for all regions
regions <- fread("qtl_project/.locusPackRat/input/regions.csv")
for (region in regions$region_id) {
  generateLocusZoomPlot_v2(
    region_id = region,
    project_dir = "qtl_project"
  )
}
```

## Working with QTL Data

```r
# Initialize with QTL regions
qtl_results <- data.frame(
  chr = c(5, 10),
  start = c(10000000, 50000000),
  end = c(15000000, 55000000),
  peak_pos = c(12500000, 52500000),
  lod = c(8.5, 6.2),
  trait = c("heart_rate", "blood_pressure")
)

initPackRat(qtl_results, mode = "region",
           species = "mouse", genome = "mm39",
           project_dir = "heart_qtl")

# Add scan data
addRatTable(scan_data, "scans", "region", project_dir = "heart_qtl")

# Add gene priorities
addRatTable(gene_scores, "priorities", "gene", project_dir = "heart_qtl")

# Generate comprehensive report
generateGeneSheet(
  format = "excel",
  include_supplementary = TRUE,
  project_dir = "heart_qtl"
)
```

## Migration from Legacy Functions

If you're using the old packet_core workflow, temporary wrappers are available:
- `generateLocusZoomPlot_legacy()` - Wraps old generateLocusZoomPlot
- `generateGeneInfoExcel_legacy()` - Wraps old generateGeneInfoExcel

See the [package documentation](https://github.com/RauLabUNC/locusPackRat) for detailed migration examples.

## Citation

Gural B, Kimball T, Luu A, Rau CD. locusPackRat: A Flexible Framework for Prioritizing Candidate Genes from GWAS and other Gene-Level Studies. *In preparation* (2025).

## Support

- Issues: https://github.com/RauLabUNC/locusPackRat/issues
- Contact: bgural@unc.edu

## License

MIT