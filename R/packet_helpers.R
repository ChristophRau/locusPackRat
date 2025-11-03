# Internal helper functions for packet generation
# These are not exported and are only used internally by generateLocusPacket()

#' Get genes in genomic region
#'
#' @param target_chr Character, chromosome
#' @param target_start Numeric, start position in bp
#' @param target_end Numeric, end position in bp
#' @param gene_data data.table with gene coordinates
#'
#' @return data.table of genes overlapping the region
#' @keywords internal
#' @noRd
.getGenesInRegion <- function(target_chr, target_start, target_end, gene_data) {
  gene_data[chr == target_chr & start_bp < target_end & end_bp > target_start, ]
}

#' Get founder SNPs in genomic region
#'
#' @param target_chr Character, chromosome
#' @param target_start Numeric, start position in bp
#' @param target_end Numeric, end position in bp
#' @param snp_data data.table with SNP coordinates
#'
#' @return data.table of SNPs in the region
#' @keywords internal
#' @noRd
.getFounderSnpsInRegion <- function(target_chr, target_start, target_end, snp_data) {
  snp_data[CHR == target_chr & mm39 < target_end & mm39 > target_start, ]
}

#' Get default heart-related phenotype terms
#'
#' @return Character string with regex pattern for cardiac terms
#' @keywords internal
#' @noRd
.getDefaultHeartPattern <- function() {
  heart_terms <- c(
    # Terms for human assoc nomenclature
    "heart", "cardi", "coronary", "atri", "ventric", "myocard",
    "hypertroph", "fibril", "valv", "rhythm", "arrhythm", "tachycard", "bradycard",
    "blood pressure", "pulse", "aort", "angina", "stroke", "hypertens", "thromb", "pulse pressure",
    "infarct", "ischem",
    # Terms for mouse assoc nomenclature
    "heart weight", "heart atrium", "heart ventricle", "heart morphology", "myocardi",
    "cardiac", "cardio", "atrioventricular", "pericardial", "ventricular septal defect",
    "atrial", "thin myocardium", "common atrioventricular valve",
    "coronary vessel", "blood vessel", "vascular", "vasculature", "aort",
    "lymphatic vessel", "yolk sac vascular", "arch artery", "pharyngeal arch artery",
    "heart failure", "circulat", "thromb", "hemorrhage", "congest", "blood circulation",
    "response of heart", "cardiac muscle", "myocardium layer", "myocardial fiber"
  )

  # Create pattern with word boundaries
  paste0("\\b(", paste(heart_terms, collapse = "|"), ")\\b|", "LV\\.")
}

#' Extract cardiac traits from human disease summary
#'
#' @param disease_summary Character, comma-separated disease names
#' @param heart_pattern Character, regex pattern for cardiac terms
#'
#' @return Character string with matched cardiac traits
#' @keywords internal
#' @noRd
.extractHumanCardiacTraits <- function(disease_summary, heart_pattern) {
  if (is.na(disease_summary) || disease_summary == "") {
    return("No cardiac traits identified")
  }

  # Split and trim
  all_traits <- trimws(unlist(strsplit(disease_summary, ",")))

  # Filter for cardiac traits
  cardiac_traits <- all_traits[grepl(heart_pattern, all_traits, ignore.case = TRUE)]

  # Return unique
  if (length(cardiac_traits) > 0) {
    return(paste(unique(cardiac_traits), collapse = ", "))
  } else {
    return("No cardiac traits identified")
  }
}

#' Extract cardiac phenotypes from mouse phenotype summary
#'
#' @param phenotype_summary Character, comma-separated phenotype names
#' @param heart_pattern Character, regex pattern for cardiac terms
#'
#' @return Character string with matched cardiac phenotypes
#' @keywords internal
#' @noRd
.extractMouseCardiacTraits <- function(phenotype_summary, heart_pattern) {
  if (is.na(phenotype_summary) || phenotype_summary == "") {
    return("No cardiac phenotypes identified")
  }

  # Split and trim
  all_phenotypes <- trimws(unlist(strsplit(phenotype_summary, ", ")))

  # Filter for cardiac phenotypes
  cardiac_phenotypes <- all_phenotypes[grepl(heart_pattern, all_phenotypes, ignore.case = TRUE)]

  # Return unique
  if (length(cardiac_phenotypes) > 0) {
    return(paste(unique(cardiac_phenotypes), collapse = ", "))
  } else {
    return("No cardiac phenotypes identified")
  }
}

#' Format genes with human-only cardiac associations
#'
#' @param genes_df data.frame with gene annotations
#' @param heart_pattern Character, regex pattern
#'
#' @return Character string formatted for README
#' @keywords internal
#' @noRd
.formatHumanOnly <- function(genes_df, heart_pattern) {
  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$`Mouse Gene Symbol`[i]
    traits <- .extractHumanCardiacTraits(genes_df$`Human Disease Summary`[i], heart_pattern)
    result[i] <- paste0(gene, ": ", traits)
  }
  return(paste(result, collapse = "\n        "))
}

#' Format genes with mouse-only cardiac associations
#'
#' @param genes_df data.frame with gene annotations
#' @param heart_pattern Character, regex pattern
#'
#' @return Character string formatted for README
#' @keywords internal
#' @noRd
.formatMouseOnly <- function(genes_df, heart_pattern) {
  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$`Mouse Gene Symbol`[i]
    phenotypes <- .extractMouseCardiacTraits(genes_df$`Mouse Phenotype Summary (MGI)`[i], heart_pattern)
    result[i] <- paste0(gene, ": ", phenotypes)
  }
  return(paste(result, collapse = "\n        "))
}

#' Format genes with both human and mouse cardiac associations
#'
#' @param genes_df data.frame with gene annotations
#' @param heart_pattern Character, regex pattern
#'
#' @return Character string formatted for README
#' @keywords internal
#' @noRd
.formatBoth <- function(genes_df, heart_pattern) {
  if (nrow(genes_df) == 0) return("None identified")

  result <- character(nrow(genes_df))
  for (i in 1:nrow(genes_df)) {
    gene <- genes_df$`Mouse Gene Symbol`[i]
    human_traits <- .extractHumanCardiacTraits(genes_df$`Human Disease Summary`[i], heart_pattern)
    mouse_phenotypes <- .extractMouseCardiacTraits(genes_df$`Mouse Phenotype Summary (MGI)`[i], heart_pattern)
    result[i] <- paste0(gene, ":\n            Human: ", human_traits, "\n            Mouse: ", mouse_phenotypes)
  }
  return(paste(result, collapse = "\n\n        "))
}

#' Generate cardiac gene summary README
#'
#' @param gene_table data.table with gene annotations
#' @param locus_cluster data.frame with locus info
#' @param locus_name Character, locus identifier
#' @param heart_pattern Character, regex pattern for cardiac terms
#' @param output_dir Character, directory to save README
#'
#' @return NULL (writes file)
#' @keywords internal
#' @noRd
.generateCardiacSummary <- function(gene_table, locus_cluster, locus_name,
                                    heart_pattern, output_dir) {

  # Filter for heart-related genes
  human_genes <- gene_table[grepl(heart_pattern, gene_table$`Human Disease Summary`, ignore.case = TRUE), ]
  mouse_genes <- gene_table[grepl(heart_pattern, gene_table$`Mouse Phenotype Summary (MGI)`, ignore.case = TRUE), ]

  # Get subsets
  mouse_only_genes <- mouse_genes[!mouse_genes$`Mouse Gene Symbol` %in% human_genes$`Mouse Gene Symbol`, ]
  mouse_only_genes <- mouse_only_genes[order(mouse_only_genes$`Mouse Gene Symbol`), ]

  human_only_genes <- human_genes[!human_genes$`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`, ]
  human_only_genes <- human_only_genes[order(human_only_genes$`Mouse Gene Symbol`), ]

  human_mouse_genes <- human_genes[human_genes$`Mouse Gene Symbol` %in% mouse_genes$`Mouse Gene Symbol`, ]
  human_mouse_genes <- human_mouse_genes[order(human_mouse_genes$`Mouse Gene Symbol`), ]

  # Summary statistics
  human_only_count <- nrow(human_only_genes)
  mouse_only_count <- nrow(mouse_only_genes)
  both_count <- nrow(human_mouse_genes)
  total_count <- human_only_count + mouse_only_count + both_count

  # Get trait info by drug
  ctrl_traits <- unique(locus_cluster[locus_cluster$drug == "Ctrl", "trait"])
  iso_traits <- unique(locus_cluster[locus_cluster$drug == "Iso", "trait"])

  # Create summary text
  summary_text <- paste(
    "# CARDIAC GENE SUMMARY\n",
    "Summary for Locus:", locus_name,
    "\nAssociated Trait(s):",
    "\n    Ctrl:", if (length(ctrl_traits) > 0) paste(ctrl_traits, collapse = ", ") else "None",
    "\n    Iso:", if (length(iso_traits) > 0) paste(iso_traits, collapse = ", ") else "None",
    "\nGenomic Region (chr", locus_cluster$chr[1], "):",
    min(locus_cluster$upper_pos_lod_drop), "-",
    max(locus_cluster$lower_pos_lod_drop), "(mm39)",
    "\n\n# Summary of known associations of genes with this loci",
    "\nTotal cardiac-related genes in region:", total_count,
    "\n  - Human cardiac traits only:", human_only_count,
    "\n  - Mouse cardiac phenotypes only:", mouse_only_count,
    "\n  - Both human and mouse:", both_count,
    "\n\n# DETAILED GENE ANNOTATIONS",
    "\n\n## 1. Genes known to be associated with human cardiac traits only (", human_only_count, "):",
    "\n        ", .formatHumanOnly(human_only_genes, heart_pattern),
    "\n\n## 2. Genes known to be associated with mouse cardiac phenotypes only (", mouse_only_count, "):",
    "\n        ", .formatMouseOnly(mouse_only_genes, heart_pattern),
    "\n\n## 3. Genes associated with both human and mouse cardiac traits (", both_count, "):",
    "\n        ", .formatBoth(human_mouse_genes, heart_pattern),
    "\n\n# NOTES",
    "\n- All 'Associated Traits (Drug)' in the original data are our own miQTL mappings",
    "\n- See cardiac_genes_summary.csv for a structured version of this data"
  )

  # Write to file
  writeLines(summary_text, file.path(output_dir, "README_summary.txt"))

  # Also create structured CSV
  cardiac_genes_csv <- rbind(
    if (nrow(human_only_genes) > 0) {
      data.frame(
        Gene_Symbol = human_only_genes$`Mouse Gene Symbol`,
        Ensembl_ID = human_only_genes$`Mouse Ensembl ID`,
        Category = "Human Only",
        Human_Cardiac_Traits = sapply(human_only_genes$`Human Disease Summary`,
                                      function(x) .extractHumanCardiacTraits(x, heart_pattern)),
        Mouse_Cardiac_Phenotypes = NA,
        Chromosome = human_only_genes$chr,
        Start_Position = human_only_genes$start_bp,
        End_Position = human_only_genes$end_bp,
        NRVM_CPM = human_only_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    },
    if (nrow(mouse_only_genes) > 0) {
      data.frame(
        Gene_Symbol = mouse_only_genes$`Mouse Gene Symbol`,
        Ensembl_ID = mouse_only_genes$`Mouse Ensembl ID`,
        Category = "Mouse Only",
        Human_Cardiac_Traits = NA,
        Mouse_Cardiac_Phenotypes = sapply(mouse_only_genes$`Mouse Phenotype Summary (MGI)`,
                                          function(x) .extractMouseCardiacTraits(x, heart_pattern)),
        Chromosome = mouse_only_genes$chr,
        Start_Position = mouse_only_genes$start_bp,
        End_Position = mouse_only_genes$end_bp,
        NRVM_CPM = mouse_only_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    },
    if (nrow(human_mouse_genes) > 0) {
      data.frame(
        Gene_Symbol = human_mouse_genes$`Mouse Gene Symbol`,
        Ensembl_ID = human_mouse_genes$`Mouse Ensembl ID`,
        Category = "Both",
        Human_Cardiac_Traits = sapply(human_mouse_genes$`Human Disease Summary`,
                                      function(x) .extractHumanCardiacTraits(x, heart_pattern)),
        Mouse_Cardiac_Phenotypes = sapply(human_mouse_genes$`Mouse Phenotype Summary (MGI)`,
                                          function(x) .extractMouseCardiacTraits(x, heart_pattern)),
        Chromosome = human_mouse_genes$chr,
        Start_Position = human_mouse_genes$start_bp,
        End_Position = human_mouse_genes$end_bp,
        NRVM_CPM = human_mouse_genes$`Avg NRVM CPM (Ctrl)`,
        stringsAsFactors = FALSE
      )
    }
  )

  # Save CSV
  if (!is.null(cardiac_genes_csv) && nrow(cardiac_genes_csv) > 0) {
    data.table::fwrite(cardiac_genes_csv, file.path(output_dir, "cardiac_genes_summary.csv"))
  }

  message("Generated README summary and cardiac genes CSV")
}
