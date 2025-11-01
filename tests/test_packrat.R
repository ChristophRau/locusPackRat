#!/usr/bin/env Rscript

# Test script for PackRat package
# This demonstrates the basic workflow of the package

# Install the package locally (if needed)
# devtools::install(".")

# Load the package
library(PackRat)

# Create sample data for testing
set.seed(123)

# 1. Sample gene data (as if from a GWAS locus)
genes_df <- data.frame(
  gene_id = paste0("ENSMUSG", sprintf("%08d", 1:10)),
  gene_symbol = c("Abcb10", "Acsl5", "Cisd2", "Pdlim5", "Manba",
                  "Gene6", "Gene7", "Gene8", "Gene9", "Gene10"),
  chr = rep(3, 10),
  start = seq(131000000, 132000000, length.out = 10),
  end = seq(131010000, 132010000, length.out = 10),
  biotype = c(rep("protein_coding", 8), "lincRNA", "pseudogene"),
  strand = sample(c("+", "-"), 10, replace = TRUE)
)

# 2. Sample orthology data
orthology_df <- data.frame(
  gene_id = genes_df$gene_id[1:8],  # Only protein-coding have orthologs
  human_gene_id = paste0("ENSG", sprintf("%08d", 1:8)),
  human_symbol = c("ABCB10", "ACSL5", "CISD2", "PDLIM5", "MANBA",
                   "GENE6", "GENE7", "GENE8"),
  identity_pct = runif(8, 85, 99)
)

# 3. Sample expression data (CPM values)
expression_df <- data.frame(
  gene_id = genes_df$gene_id,
  CPM_heart = c(150.2, 89.5, 234.1, 178.9, 45.3,
                12.4, 189.2, 67.8, 3.2, 0.5),
  CPM_liver = c(45.3, 234.5, 123.4, 34.2, 189.3,
                234.1, 45.2, 123.4, 5.6, 1.2)
)

# 4. Sample eQTL data
eqtl_df <- data.frame(
  gene_id = genes_df$gene_id[c(3, 4, 7)],  # Only some genes have eQTLs
  snp_id = c("rs123456", "rs234567", "rs345678"),
  p_value = c(1e-8, 5e-6, 3e-5),
  beta = c(0.45, -0.32, 0.28)
)

# 5. Sample variant data
variants_df <- data.frame(
  gene_id = rep(genes_df$gene_id[c(1, 3, 4, 5)], each = 2),
  variant_type = c("missense", "synonymous", "missense", "stop_gained",
                   "missense", "missense", "frameshift", "missense"),
  position = c(131000100, 131000200, 131200100, 131200200,
               131300100, 131300200, 131400100, 131400200),
  strain = rep(c("PWK/PhJ", "CAST/EiJ"), 4),
  impact = c("moderate", "low", "moderate", "high",
             "moderate", "moderate", "high", "moderate")
)

# 6. Sample phenotype data
phenotype_df <- data.frame(
  gene_id = rep(genes_df$gene_id[1:5], each = 2),
  phenotype = c("abnormal cardiac muscle contractility",
                "increased heart weight",
                "abnormal skeletal muscle morphology",
                "decreased cardiac muscle contractility",
                "dilated cardiomyopathy",
                "abnormal heart morphology",
                "abnormal enzyme activity",
                "decreased circulating glucose level",
                "abnormal cardiac ventricle morphology",
                "increased susceptibility to cardiac hypertrophy"),
  source = rep(c("MGI", "IMPC"), 5),
  p_value = runif(10, 1e-10, 1e-3)
)

cat("Sample data created successfully\\n\\n")

# Test 1: Build comprehensive gene table
cat("Test 1: Building comprehensive gene table...\\n")
gene_table <- build_gene_table(
  genes_df,
  orthology_df = orthology_df,
  expression_df = expression_df,
  eqtl_df = eqtl_df,
  variants_df = variants_df
)
print(str(gene_table))
cat("\\n")

# Test 2: Filter genes
cat("Test 2: Filtering genes...\\n")
filtered_genes <- filter_genes(
  gene_table,
  require_human_ortholog = TRUE,
  min_expression = 50,
  biotype_filter = "protein_coding"
)
cat("Filtered from", nrow(gene_table), "to", nrow(filtered_genes), "genes\\n")
cat("Filtered genes:", paste(filtered_genes$gene_symbol, collapse = ", "), "\\n\\n")

# Test 3: Extract phenotypes
cat("Test 3: Extracting cardiac-related phenotypes...\\n")
cardiac_phenotypes <- extract_phenotypes(
  phenotype_df,
  use_sets = c("cardiac")
)
cat("Found", nrow(cardiac_phenotypes), "cardiac phenotypes\\n\\n")

# Test 4: Score phenotype relevance
cat("Test 4: Scoring phenotype relevance...\\n")
scored_phenotypes <- score_phenotype_relevance(
  cardiac_phenotypes,
  primary_keywords = c("cardiac", "heart"),
  secondary_keywords = c("muscle", "contractility")
)
cat("Top phenotype:", scored_phenotypes$phenotype[1],
    "with score:", scored_phenotypes$relevance_score[1], "\\n\\n")

# Test 5: Create Excel workbook
cat("Test 5: Creating Excel workbook...\\n")
output_file <- "test_locus_chr3.xlsx"
create_gene_workbook(
  gene_table,
  output_file = output_file
)
cat("Excel workbook created:", output_file, "\\n\\n")

# Test 6: Generate locus report
cat("Test 6: Generating locus report...\\n")
locus_info <- list(
  chr = 3,
  start = 131000000,
  end = 132000000,
  trait = "Cardiac hypertrophy response",
  p_value = 1e-12,
  n_genes = nrow(genes_df)
)

report <- generate_locus_report(
  locus_info,
  gene_table,
  output_format = "text"
)
cat(substr(report, 1, 500), "...\\n\\n")  # Show first 500 chars

# Test 7: Identify top candidates
cat("Test 7: Identifying top candidate genes...\\n")
top_candidates <- identify_top_candidates(
  gene_table,
  prioritization_criteria = list(
    require_eqtl = FALSE,
    require_variant = TRUE,
    min_expression = 50
  ),
  n_top = 3
)
cat("Top candidates:\\n")
for(i in 1:nrow(top_candidates)) {
  cat("  ", i, ". ", top_candidates$gene_symbol[i],
      " (score: ", round(top_candidates$priority_score[i], 2), ")\\n", sep = "")
}

cat("\\n✓ All tests completed successfully!\\n")