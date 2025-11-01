#!/usr/bin/env Rscript

#' Demonstration of relational filtering in PackRat
#'
#' This shows how filterGenes works as a general relational filtering tool

library(PackRat)
set.seed(123)

cat("=====================================================\n")
cat("RELATIONAL FILTERING DEMONSTRATION\n")
cat("=====================================================\n\n")

# Create sample data tables
cat("Creating sample data tables...\n\n")

# Main gene table
genes <- data.frame(
  gene_id = paste0("ENSMUSG", sprintf("%08d", 1:20)),
  symbol = paste0("Gene", 1:20),
  chr = sample(1:5, 20, replace = TRUE),
  biotype = sample(c("protein_coding", "lincRNA", "pseudogene"), 20,
                   replace = TRUE, prob = c(0.7, 0.2, 0.1)),
  stringsAsFactors = FALSE
)

# Expression data (separate table)
expression <- data.frame(
  gene_id = genes$gene_id,
  TPM_heart = runif(20, 0, 500),
  TPM_liver = runif(20, 0, 300),
  TPM_brain = runif(20, 0, 400),
  stringsAsFactors = FALSE
)

# GWAS hits (subset of genes)
gwasHits <- data.frame(
  gene_id = sample(genes$gene_id, 8),
  pvalue = 10^-runif(8, 5, 12),
  trait = sample(c("cardiac", "metabolic"), 8, replace = TRUE),
  stringsAsFactors = FALSE
)

# Ortholog mapping (different column names)
orthologs <- data.frame(
  mouse_gene = genes$gene_id[1:15],  # Not all genes have orthologs
  human_gene = paste0("ENSG", sprintf("%08d", 1:15)),
  identity = runif(15, 70, 99),
  stringsAsFactors = FALSE
)

# Disease associations (by symbol, not ID)
diseaseGenes <- data.frame(
  gene_symbol = sample(genes$symbol, 10),
  disease = sample(c("cardiomyopathy", "diabetes", "obesity"), 10, replace = TRUE),
  evidence_score = runif(10, 0.1, 0.9),
  stringsAsFactors = FALSE
)

cat("Created tables:\n")
cat("  - genes:", nrow(genes), "genes\n")
cat("  - expression:", nrow(expression), "measurements\n")
cat("  - gwasHits:", nrow(gwasHits), "significant hits\n")
cat("  - orthologs:", nrow(orthologs), "mouse-human pairs\n")
cat("  - diseaseGenes:", nrow(diseaseGenes), "disease associations\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 1: Simple semi-join (keep genes in a list)
# -----------------------------------------------------------------------------
cat("========================================\n")
cat("EXAMPLE 1: Filter by GWAS hits\n")
cat("========================================\n")
cat("Goal: Keep only genes that are GWAS hits\n\n")

result1 <- filterGenes(
  inputTable = genes,
  referenceTable = gwasHits,
  by = "gene_id",
  joinType = "semi",
  verbose = TRUE
)

cat("\nGenes in GWAS hits:", paste(result1$symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 2: Anti-join (exclude genes)
# -----------------------------------------------------------------------------
cat("========================================\n")
cat("EXAMPLE 2: Exclude disease genes\n")
cat("========================================\n")
cat("Goal: Remove genes associated with diseases\n\n")

result2 <- filterGenes(
  inputTable = genes,
  referenceTable = diseaseGenes,
  by = c("symbol" = "gene_symbol"),  # Different column names
  joinType = "anti",
  verbose = TRUE
)

cat("\nGenes NOT in disease list:", paste(result2$symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 3: Join and filter
# -----------------------------------------------------------------------------
cat("========================================\n")
cat("EXAMPLE 3: Join with expression and filter\n")
cat("========================================\n")
cat("Goal: Add expression data and filter for highly expressed genes\n\n")

result3 <- filterGenes(
  inputTable = genes,
  referenceTable = expression,
  by = "gene_id",
  joinType = "left",
  filters = list(
    makeFilter("TPM_heart", ">", 250),
    makeFilter("biotype", "==", "protein_coding")
  ),
  verbose = TRUE
)

cat("\nHighly expressed protein-coding genes:\n")
for (i in 1:nrow(result3)) {
  cat(sprintf("  %s: %.1f TPM\n", result3$symbol[i], result3$TPM_heart[i]))
}

# -----------------------------------------------------------------------------
# EXAMPLE 4: Multi-step relational filtering
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("EXAMPLE 4: Complex multi-table workflow\n")
cat("========================================\n")
cat("Goal: Find GWAS hits with human orthologs and high expression\n\n")

# Step 1: Get GWAS hits with expression data
cat("Step 1: Join GWAS hits with expression\n")
step1 <- filterGenes(
  inputTable = gwasHits,
  referenceTable = expression,
  by = "gene_id",
  joinType = "left",
  verbose = FALSE
)
cat("  ->", nrow(step1), "GWAS hits with expression data\n")

# Step 2: Filter for cardiac GWAS hits with high expression
cat("Step 2: Filter for cardiac hits with high heart expression\n")
step2 <- filterGenes(
  inputTable = step1,
  filters = list(
    makeFilter("trait", "==", "cardiac"),
    makeFilter("TPM_heart", ">", 100)
  ),
  verbose = FALSE
)
cat("  ->", nrow(step2), "cardiac GWAS hits with high expression\n")

# Step 3: Check for human orthologs
cat("Step 3: Keep only those with human orthologs\n")
step3 <- filterGenes(
  inputTable = step2,
  referenceTable = orthologs,
  by = c("gene_id" = "mouse_gene"),
  joinType = "inner",
  verbose = FALSE
)
cat("  ->", nrow(step3), "genes meet all criteria\n")

cat("\nFinal candidates:\n")
if (nrow(step3) > 0) {
  for (i in 1:nrow(step3)) {
    cat(sprintf("  %s -> %s (%.1f%% identity, %.1f TPM, p=%.2e)\n",
                step3$gene_id[i], step3$human_gene[i],
                step3$identity[i], step3$TPM_heart[i], step3$pvalue[i]))
  }
} else {
  cat("  No genes met all criteria\n")
}

# -----------------------------------------------------------------------------
# EXAMPLE 5: Using custom functions with relational filtering
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("EXAMPLE 5: Combine relational and custom filters\n")
cat("========================================\n")
cat("Goal: GWAS hits with tissue-specific expression\n\n")

result5 <- filterGenes(
  inputTable = genes,
  referenceTable = gwasHits,
  by = "gene_id",
  joinType = "inner",  # Get GWAS hits with all their info
  filters = list(
    # Custom function for tissue specificity
    function(df) {
      # Join with expression data first
      df_expr <- merge(df, expression, by = "gene_id")
      # Calculate tissue specificity (highest / mean of others)
      heart_specific <- df_expr$TPM_heart > 2 * mean(c(df_expr$TPM_liver, df_expr$TPM_brain))
      return(heart_specific)
    }
  ),
  verbose = TRUE
)

cat("\nHeart-specific GWAS hits:", paste(result5$symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 6: Building gene sets
# -----------------------------------------------------------------------------
cat("========================================\n")
cat("EXAMPLE 6: Building custom gene sets\n")
cat("========================================\n\n")

# Create multiple gene sets
highExpression <- filterGenes(
  inputTable = expression,
  filters = list(
    function(df) rowSums(df[, c("TPM_heart", "TPM_liver", "TPM_brain")] > 200) >= 2
  ),
  keepColumns = "gene_id",
  verbose = FALSE
)

hasOrtholog <- filterGenes(
  inputTable = orthologs,
  filters = list(makeFilter("identity", ">", 80)),
  keepColumns = "mouse_gene",
  verbose = FALSE
)

# Find intersection of gene sets
cat("High expression genes:", nrow(highExpression), "\n")
cat("Genes with good orthologs:", nrow(hasOrtholog), "\n")

intersection <- filterGenes(
  inputTable = highExpression,
  referenceTable = hasOrtholog,
  by = c("gene_id" = "mouse_gene"),
  joinType = "semi",
  verbose = FALSE
)

cat("Intersection:", nrow(intersection), "genes\n\n")

# -----------------------------------------------------------------------------
cat("=====================================================\n")
cat("KEY INSIGHTS\n")
cat("=====================================================\n\n")

cat("1. filterGenes is a general relational filtering tool\n")
cat("   - Works with ANY data frames, not just genes\n")
cat("   - Handles different column names via 'by' parameter\n\n")

cat("2. Join types provide different filtering strategies:\n")
cat("   - semi: Keep rows that match (like %in%)\n")
cat("   - anti: Keep rows that DON'T match (like !%in%)\n")
cat("   - inner: Combine matching rows from both tables\n")
cat("   - left: Keep all input rows, add reference columns\n\n")

cat("3. Combine relational and value-based filtering:\n")
cat("   - First join tables to establish relationships\n")
cat("   - Then filter on values from either table\n\n")

cat("4. Chain operations for complex workflows:\n")
cat("   - Each filterGenes call returns a data frame\n")
cat("   - Results can be piped into the next step\n")
cat("   - Build complex queries from simple operations\n\n")

cat("=====================================================\n")
cat("DEMO COMPLETE\n")
cat("=====================================================\n")