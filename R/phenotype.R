#' Extract phenotypes by keyword matching
#' 
#' @description
#' Flexible phenotype extraction from various databases (MGI, IMPC, Open Targets)
#' using keyword matching. Supports custom keyword sets for different research areas.
#' 
#' @param phenotype_df Data frame with phenotype descriptions
#' @param keywords Character vector of keywords to search for
#' @param keyword_sets Named list of predefined keyword sets (e.g., "cardiac", "metabolic")
#' @param use_sets Character vector of keyword set names to use
#' @param text_columns Columns to search for keywords (default: auto-detect)
#' @param case_sensitive Logical, should matching be case sensitive (default: FALSE)
#' @param match_type How to match keywords: "any", "all", "exact" (default: "any")
#' 
#' @return Data frame with matched phenotypes and match details
#' 
#' @examples
#' # Using custom keywords
#' pheno_df <- data.frame(
#'   gene = c("Abcb10", "Acsl5"),
#'   phenotype = c("cardiac hypertrophy", "insulin resistance")
#' )
#' cardiac_phenos <- extract_phenotypes(pheno_df, keywords = c("cardiac", "heart"))
#' 
#' # Using predefined keyword sets
#' results <- extract_phenotypes(pheno_df, use_sets = c("cardiac", "metabolic"))
#' 
#' @export
extract_phenotypes <- function(phenotype_df,
                              keywords = NULL,
                              keyword_sets = NULL,
                              use_sets = NULL,
                              text_columns = NULL,
                              case_sensitive = FALSE,
                              match_type = "any") {
  
  # Default keyword sets if none provided
  if (is.null(keyword_sets)) {
    keyword_sets <- list(
      cardiac = c("cardiac", "heart", "cardiom", "atri", "ventric", 
                 "arrhythm", "contracti", "ejection", "diastol", "systol",
                 "myocard", "pericard", "endocard", "valve", "coronary"),
      
      metabolic = c("glucose", "insulin", "diabet", "lipid", "cholesterol",
                   "triglyceride", "adipos", "obesity", "metabol", "glycem",
                   "leptin", "adiponectin", "fatty acid"),
      
      neurological = c("neuro", "brain", "cereb", "cognit", "memory", 
                      "learning", "behavior", "anxiety", "depression",
                      "motor", "sensory", "seizure", "epilep"),
      
      immunological = c("immun", "inflamm", "cytokine", "lymph", "macrophage",
                       "neutrophil", "antibody", "antigen", "allerg",
                       "autoimmun", "infection"),
      
      renal = c("kidney", "renal", "nephro", "glomerul", "tubul",
               "proteinuria", "albumin", "creatin", "urea"),
      
      hepatic = c("liver", "hepat", "cirrhosis", "steatos", "fibros",
                 "biliary", "bile", "bilirubin", "transaminase"),
      
      respiratory = c("lung", "pulmon", "respirat", "bronch", "alveol",
                     "asthma", "fibrosis", "emphysema", "pneumo"),
      
      skeletal = c("bone", "skeletal", "oste", "cartilage", "joint",
                  "arthr", "fracture", "mineral density", "growth plate"),
      
      reproductive = c("reproduct", "fertil", "gonad", "testis", "ovary",
                      "sperm", "oocyte", "pregnan", "embryo", "mating")
    )
  }
  
  # Compile keywords from sets and individual keywords
  all_keywords <- character()
  
  if (!is.null(use_sets)) {
    for (set_name in use_sets) {
      if (set_name %in% names(keyword_sets)) {
        all_keywords <- c(all_keywords, keyword_sets[[set_name]])
      } else {
        warning(sprintf("Keyword set '%s' not found", set_name))
      }
    }
  }
  
  if (!is.null(keywords)) {
    all_keywords <- c(all_keywords, keywords)
  }
  
  if (length(all_keywords) == 0) {
    stop("No keywords provided. Specify 'keywords' or 'use_sets'")
  }
  
  # Remove duplicates
  all_keywords <- unique(all_keywords)
  
  # Identify text columns to search
  if (is.null(text_columns)) {
    # Auto-detect text columns
    text_columns <- names(phenotype_df)[sapply(phenotype_df, is.character)]
    # Prioritize likely phenotype columns
    pheno_cols <- text_columns[grepl("phenotype|description|term|name", 
                                    text_columns, ignore.case = TRUE)]
    if (length(pheno_cols) > 0) {
      text_columns <- pheno_cols
    }
  }
  
  if (length(text_columns) == 0) {
    stop("No text columns found to search")
  }
  
  # Create combined text field for searching
  phenotype_df$search_text <- apply(phenotype_df[text_columns], 1, 
                                   function(x) paste(x, collapse = " "))
  
  # Apply case sensitivity
  if (!case_sensitive) {
    phenotype_df$search_text <- tolower(phenotype_df$search_text)
    all_keywords <- tolower(all_keywords)
  }
  
  # Perform matching based on match_type
  if (match_type == "any") {
    # Match any keyword
    pattern <- paste0("(", paste(all_keywords, collapse = "|"), ")")
    phenotype_df$matched <- grepl(pattern, phenotype_df$search_text)
    
    # Find which keywords matched
    phenotype_df$matched_keywords <- NA
    for (i in which(phenotype_df$matched)) {
      matched_kw <- all_keywords[sapply(all_keywords, 
                                       function(kw) grepl(kw, phenotype_df$search_text[i]))]
      phenotype_df$matched_keywords[i] <- paste(matched_kw, collapse = "; ")
    }
    
  } else if (match_type == "all") {
    # Must match all keywords
    phenotype_df$matched <- apply(sapply(all_keywords, 
                                        function(kw) grepl(kw, phenotype_df$search_text)), 
                                 1, all)
    phenotype_df$matched_keywords <- ifelse(phenotype_df$matched, 
                                           paste(all_keywords, collapse = "; "), 
                                           NA)
    
  } else if (match_type == "exact") {
    # Exact word matching
    phenotype_df$matched <- FALSE
    phenotype_df$matched_keywords <- NA
    
    for (kw in all_keywords) {
      pattern <- paste0("\\b", kw, "\\b")
      matches <- grepl(pattern, phenotype_df$search_text)
      phenotype_df$matched <- phenotype_df$matched | matches
      
      for (i in which(matches)) {
        if (is.na(phenotype_df$matched_keywords[i])) {
          phenotype_df$matched_keywords[i] <- kw
        } else {
          phenotype_df$matched_keywords[i] <- paste(phenotype_df$matched_keywords[i], 
                                                   kw, sep = "; ")
        }
      }
    }
  }
  
  # Remove temporary search_text column
  phenotype_df$search_text <- NULL
  
  # Add match score (number of keywords matched)
  phenotype_df$match_score <- sapply(strsplit(phenotype_df$matched_keywords, "; "), 
                                    function(x) length(x[!is.na(x)]))
  phenotype_df$match_score[is.na(phenotype_df$matched_keywords)] <- 0
  
  # Filter to matched records
  result <- phenotype_df[phenotype_df$matched, ]
  result$matched <- NULL  # Remove redundant column
  
  # Sort by match score
  result <- result[order(result$match_score, decreasing = TRUE), ]
  
  return(result)
}

#' Score phenotype relevance
#' 
#' @description
#' Calculate relevance scores for phenotypes based on multiple criteria.
#' 
#' @param phenotype_df Data frame with phenotype matches
#' @param score_weights Named list of weights for different scoring components
#' @param normalize Normalize scores to 0-1 range (default: TRUE)
#' 
#' @return Data frame with relevance scores added
#' 
#' @export
score_phenotype_relevance <- function(phenotype_df,
                                     score_weights = NULL,
                                     normalize = TRUE) {
  
  # Default weights if not provided
  if (is.null(score_weights)) {
    score_weights <- list(
      keyword_matches = 1.0,    # Number of keywords matched
      database_score = 0.5,      # Database-provided confidence score
      evidence_count = 0.3,      # Number of supporting studies
      species_match = 0.8        # Same species as study
    )
  }
  
  phenotype_df$relevance_score <- 0
  
  # Add component scores
  if ("match_score" %in% names(phenotype_df) && "keyword_matches" %in% names(score_weights)) {
    phenotype_df$relevance_score <- phenotype_df$relevance_score + 
      phenotype_df$match_score * score_weights$keyword_matches
  }
  
  if ("database_score" %in% names(phenotype_df) && "database_score" %in% names(score_weights)) {
    phenotype_df$relevance_score <- phenotype_df$relevance_score + 
      phenotype_df$database_score * score_weights$database_score
  }
  
  if ("evidence_count" %in% names(phenotype_df) && "evidence_count" %in% names(score_weights)) {
    phenotype_df$relevance_score <- phenotype_df$relevance_score + 
      phenotype_df$evidence_count * score_weights$evidence_count
  }
  
  # Normalize if requested
  if (normalize && max(phenotype_df$relevance_score) > 0) {
    phenotype_df$relevance_score <- phenotype_df$relevance_score / 
      max(phenotype_df$relevance_score)
  }
  
  return(phenotype_df)
}

#' Summarize phenotypes by gene
#' 
#' @description
#' Aggregate phenotype information at the gene level.
#' 
#' @param phenotype_df Data frame with phenotype information
#' @param gene_col Column name for gene identifier
#' @param collapse_sep Separator for collapsed text fields
#' 
#' @return Summary data frame with one row per gene
#' 
#' @export
summarize_phenotypes_by_gene <- function(phenotype_df,
                                        gene_col = "gene",
                                        collapse_sep = "; ") {
  
  if (!gene_col %in% names(phenotype_df)) {
    stop(sprintf("Gene column '%s' not found", gene_col))
  }
  
  # Group by gene and summarize
  summary_df <- aggregate(
    phenotype_df,
    by = list(gene = phenotype_df[[gene_col]]),
    FUN = function(x) {
      if (is.numeric(x)) {
        # For numeric columns, take max or mean
        if (all(x == round(x))) {
          return(sum(x))  # Counts
        } else {
          return(mean(x))  # Scores
        }
      } else {
        # For text, collapse unique values
        unique_vals <- unique(x[!is.na(x)])
        if (length(unique_vals) > 3) {
          return(paste(c(unique_vals[1:3], "..."), collapse = collapse_sep))
        } else {
          return(paste(unique_vals, collapse = collapse_sep))
        }
      }
    }
  )
  
  # Add count of phenotypes per gene
  pheno_counts <- table(phenotype_df[[gene_col]])
  summary_df$n_phenotypes <- pheno_counts[match(summary_df$gene, names(pheno_counts))]
  
  return(summary_df)
}