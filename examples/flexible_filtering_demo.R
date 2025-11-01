#!/usr/bin/env Rscript

#' Demonstration of flexible filtering approach for PackRat
#'
#' This script shows how the new flexible filtering system works
#' compared to the old hardcoded approach

library(PackRat)

# Create sample gene data
genes_df <- data.frame(
  gene_id = paste0("ENSMUSG", sprintf("%08d", 1:20)),
  gene_symbol = paste0("Gene", 1:20),
  chr = c(rep(3, 10), rep(5, 10)),
  biotype = c(rep("protein_coding", 15), rep("lincRNA", 3), rep("pseudogene", 2)),
  TPM_heart = runif(20, 0, 500),
  TPM_liver = runif(20, 0, 300),
  TPM_brain = runif(20, 0, 400),
  has_human_ortholog = c(rep(TRUE, 15), rep(FALSE, 5)),
  ortholog_identity = c(runif(15, 70, 99), rep(NA, 5)),
  has_eqtl = sample(c(TRUE, FALSE), 20, replace = TRUE, prob = c(0.3, 0.7)),
  eqtl_pvalue = ifelse(genes_df$has_eqtl, 10^-runif(20, 3, 10), NA),
  n_variants = sample(0:5, 20, replace = TRUE),
  has_deleterious = genes_df$n_variants > 0 & sample(c(TRUE, FALSE), 20, replace = TRUE),
  differential_expr = sample(c(TRUE, FALSE), 20, replace = TRUE, prob = c(0.4, 0.6)),
  log2FC = ifelse(genes_df$differential_expr, rnorm(20, 0, 2), 0),
  stringsAsFactors = FALSE
)

cat("===========================================\n")
cat("FLEXIBLE FILTERING DEMONSTRATION FOR PACKRAT\n")
cat("===========================================\n\n")

cat("Sample data created with", nrow(genes_df), "genes\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 1: Basic filtering
# -----------------------------------------------------------------------------
cat("EXAMPLE 1: Basic Filtering\n")
cat("--------------------------\n")
cat("Goal: Find protein-coding genes with high heart expression\n\n")

result1 <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("biotype", "==", "protein_coding"),
    make_filter("TPM_heart", ">", 200)
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result1$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 2: Multiple tissue expression
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 2: Multi-tissue Expression\n")
cat("-----------------------------------\n")
cat("Goal: Find genes expressed in heart OR liver (>100 TPM)\n\n")

result2 <- filter_genes_flexible(
  genes_df,
  filters = list(
    # Custom function for OR logic within a filter
    function(df) df$TPM_heart > 100 | df$TPM_liver > 100
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result2$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 3: Complex criteria with genetic evidence
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 3: Genes with Genetic Evidence\n")
cat("----------------------------------------\n")
cat("Goal: Protein-coding genes with eQTL or deleterious variants\n\n")

result3 <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("biotype", "==", "protein_coding"),
    function(df) df$has_eqtl | df$has_deleterious
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result3$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 4: Using presets
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 4: Using Preset Filters\n")
cat("---------------------------------\n")
cat("Goal: Use 'validated' preset plus custom chromosome filter\n\n")

# First need to adjust column names to match preset expectations
genes_df_adjusted <- genes_df
names(genes_df_adjusted)[names(genes_df_adjusted) == "TPM_heart"] <- "CPM_heart"
names(genes_df_adjusted)[names(genes_df_adjusted) == "has_human_ortholog"] <- "human_gene_id"
genes_df_adjusted$human_gene_id[genes_df_adjusted$human_gene_id == FALSE] <- NA
genes_df_adjusted$human_gene_id[genes_df_adjusted$human_gene_id == TRUE] <- paste0("ENSG", 1:sum(genes_df_adjusted$human_gene_id == TRUE, na.rm = TRUE))

result4 <- filter_genes_preset(
  genes_df_adjusted,
  preset = "validated",
  custom_filters = list(
    make_filter("chr", "==", 3)
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result4$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 5: Pattern matching
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 5: Pattern Matching\n")
cat("-----------------------------\n")
cat("Goal: Find genes with symbols starting with 'Gene1'\n\n")

result5 <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("gene_symbol", "matches", "^Gene1")
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result5$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 6: Range filtering
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 6: Range Filtering\n")
cat("---------------------------\n")
cat("Goal: Genes with moderate ortholog identity (80-95%)\n\n")

result6 <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("ortholog_identity", "between", c(80, 95))
  ),
  verbose = TRUE
)

cat("\nGenes found:", paste(result6$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 7: Combining filter modes
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 7: ANY vs ALL Filter Modes\n")
cat("------------------------------------\n")

cat("\n7a. ALL mode (default) - must meet all criteria:\n")
result7a <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("biotype", "==", "protein_coding"),
    make_filter("has_eqtl", "==", TRUE),
    make_filter("n_variants", ">", 0)
  ),
  filter_mode = "all",
  verbose = TRUE
)
cat("Genes found:", paste(result7a$gene_symbol, collapse = ", "), "\n")

cat("\n7b. ANY mode - must meet at least one criterion:\n")
result7b <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("has_eqtl", "==", TRUE),
    make_filter("n_variants", ">", 2),
    make_filter("differential_expr", "==", TRUE)
  ),
  filter_mode = "any",
  verbose = TRUE
)
cat("Genes found:", paste(result7b$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# EXAMPLE 8: Real-world prioritization workflow
# -----------------------------------------------------------------------------
cat("\nEXAMPLE 8: Complete Prioritization Workflow\n")
cat("--------------------------------------------\n")
cat("Goal: Multi-step filtering for candidate gene prioritization\n\n")

# Step 1: Initial broad filter
cat("Step 1: Get all expressed protein-coding genes\n")
step1 <- filter_genes_flexible(
  genes_df,
  filters = list(
    make_filter("biotype", "==", "protein_coding"),
    make_filter("TPM_heart", ">", 10)  # Minimal expression
  ),
  verbose = FALSE
)
cat("  ->", nrow(step1), "genes remain\n")

# Step 2: Require human ortholog for translatability
cat("Step 2: Filter for human orthologs\n")
step2 <- filter_genes_flexible(
  step1,
  filters = list(
    make_filter("has_human_ortholog", "==", TRUE),
    make_filter("ortholog_identity", ">", 75)
  ),
  verbose = FALSE
)
cat("  ->", nrow(step2), "genes remain\n")

# Step 3: Require some genetic evidence
cat("Step 3: Require genetic or expression evidence\n")
step3 <- filter_genes_flexible(
  step2,
  filters = list(
    function(df) df$has_eqtl | df$n_variants > 0 | df$differential_expr
  ),
  verbose = FALSE
)
cat("  ->", nrow(step3), "genes remain\n")

# Step 4: Prioritize by multiple evidence types
cat("Step 4: High-priority candidates (stringent)\n")
step4 <- filter_genes_flexible(
  step3,
  filters = list(
    function(df) {
      evidence_count <- as.numeric(df$has_eqtl) +
                       as.numeric(df$n_variants > 0) +
                       as.numeric(df$differential_expr)
      evidence_count >= 2  # At least 2 types of evidence
    }
  ),
  verbose = FALSE
)
cat("  ->", nrow(step4), "genes remain\n")

cat("\nFinal candidates:", paste(step4$gene_symbol, collapse = ", "), "\n\n")

# -----------------------------------------------------------------------------
# COMPARISON WITH OLD APPROACH
# -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("COMPARISON: Flexible vs Hardcoded Approach\n")
cat("===========================================\n\n")

cat("OLD APPROACH (hardcoded):\n")
cat("-----------------------\n")
cat("filter_genes(data,\n")
cat("  require_human_ortholog = TRUE,\n")
cat("  min_expression = 50,\n")
cat("  require_protein_coding = TRUE)\n")
cat("\nProblems:\n")
cat("- Fixed column names (what if your data uses TPM not CPM?)\n")
cat("- Limited filter options\n")
cat("- Can't combine conditions flexibly\n")
cat("- No support for multiple tissues\n\n")

cat("NEW APPROACH (flexible):\n")
cat("-----------------------\n")
cat("filter_genes_flexible(data,\n")
cat("  filters = list(\n")
cat("    make_filter('biotype', '==', 'protein_coding'),\n")
cat("    make_filter('TPM_heart', '>', 50),\n")
cat("    make_filter('has_human_ortholog', '==', TRUE),\n")
cat("    # Can add ANY additional criteria\n")
cat("    make_filter('ortholog_identity', '>', 85),\n")
cat("    function(df) df$TPM_heart > 100 | df$TPM_liver > 100\n")
cat("  ))\n")
cat("\nAdvantages:\n")
cat("- Works with any column names\n")
cat("- Unlimited filter combinations\n")
cat("- Custom logic with functions\n")
cat("- Can be configured from files\n")
cat("- Reusable filter specifications\n\n")

cat("===========================================\n")
cat("DEMO COMPLETE\n")
cat("===========================================\n")