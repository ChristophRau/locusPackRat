#' Flexible gene filtering with dynamic criteria
#'
#' @description
#' Filter genes using flexible, user-defined criteria on any columns.
#' Supports multiple filter types including numeric comparisons,
#' value matching, pattern matching, and custom functions.
#'
#' @param data Data frame or data.table to filter
#' @param filters List of filter specifications. Each filter can be:
#'   - A list with 'column', 'condition', and 'value' elements
#'   - A function that takes the data and returns a logical vector
#'   - A character string to be evaluated as an expression
#' @param filter_mode How to combine multiple filters: "all" (AND) or "any" (OR)
#' @param keep_columns Optional character vector of columns to retain in output
#' @param verbose Print filtering progress and statistics
#'
#' @return Filtered data frame/data.table
#'
#' @examples
#' # Example 1: Simple numeric filtering
#' filtered <- filter_genes_flexible(
#'   gene_table,
#'   filters = list(
#'     list(column = "CPM_heart", condition = ">", value = 50),
#'     list(column = "biotype", condition = "==", value = "protein_coding")
#'   )
#' )
#'
#' # Example 2: Using %in% for multiple values
#' filtered <- filter_genes_flexible(
#'   gene_table,
#'   filters = list(
#'     list(column = "biotype", condition = "%in%",
#'          value = c("protein_coding", "lincRNA")),
#'     list(column = "has_eqtl", condition = "==", value = TRUE)
#'   )
#' )
#'
#' # Example 3: Pattern matching with regex
#' filtered <- filter_genes_flexible(
#'   gene_table,
#'   filters = list(
#'     list(column = "gene_symbol", condition = "matches", value = "^Abc"),
#'     list(column = "strand", condition = "==", value = "+")
#'   )
#' )
#'
#' # Example 4: Complex filtering with custom function
#' filtered <- filter_genes_flexible(
#'   gene_table,
#'   filters = list(
#'     # Standard filters
#'     list(column = "biotype", condition = "==", value = "protein_coding"),
#'     # Custom function for complex logic
#'     function(df) df$CPM_heart > 100 | df$CPM_liver > 100,
#'     # String expression
#'     "!is.na(human_gene_id)"
#'   )
#' )
#'
#' # Example 5: Filter with NA handling
#' filtered <- filter_genes_flexible(
#'   gene_table,
#'   filters = list(
#'     list(column = "human_gene_id", condition = "not_na"),
#'     list(column = "n_variants", condition = ">=", value = 1, na.rm = TRUE)
#'   )
#' )
#'
#' @export
filter_genes_flexible <- function(data,
                                 filters = list(),
                                 filter_mode = "all",
                                 keep_columns = NULL,
                                 verbose = FALSE) {

  # Convert to data.table if needed
  if (!inherits(data, "data.table")) {
    data <- data.table::as.data.table(data)
  } else {
    data <- data.table::copy(data)
  }

  if (length(filters) == 0) {
    if (verbose) cat("No filters specified, returning all data\n")
    return(data)
  }

  # Store original row count
  original_rows <- nrow(data)

  # Initialize filter results
  filter_results <- list()

  # Process each filter
  for (i in seq_along(filters)) {
    filter <- filters[[i]]

    if (is.function(filter)) {
      # Custom function filter
      if (verbose) cat(sprintf("Applying custom function filter %d\n", i))
      filter_results[[i]] <- filter(data)

    } else if (is.character(filter) && length(filter) == 1) {
      # Expression string filter
      if (verbose) cat(sprintf("Evaluating expression: %s\n", filter))
      filter_results[[i]] <- eval(parse(text = filter), envir = data)

    } else if (is.list(filter)) {
      # Structured filter
      column <- filter$column
      condition <- filter$condition
      value <- filter$value
      na.rm <- filter$na.rm %||% FALSE

      if (!column %in% names(data)) {
        warning(sprintf("Column '%s' not found in data, skipping filter", column))
        filter_results[[i]] <- rep(TRUE, nrow(data))
        next
      }

      if (verbose) {
        cat(sprintf("Filtering %s %s %s\n", column, condition,
                   paste(value, collapse = ", ")))
      }

      # Apply condition
      result <- switch(condition,
        # Equality and inequality
        "==" = data[[column]] == value,
        "!=" = data[[column]] != value,

        # Numeric comparisons
        ">" = data[[column]] > value,
        ">=" = data[[column]] >= value,
        "<" = data[[column]] < value,
        "<=" = data[[column]] <= value,

        # Set membership
        "%in%" = data[[column]] %in% value,
        "%nin%" = !(data[[column]] %in% value),

        # Pattern matching
        "matches" = grepl(value, data[[column]], perl = TRUE),
        "not_matches" = !grepl(value, data[[column]], perl = TRUE),

        # NA checks
        "is_na" = is.na(data[[column]]),
        "not_na" = !is.na(data[[column]]),

        # Range checks
        "between" = {
          if (length(value) != 2) stop("'between' requires two values")
          data[[column]] >= value[1] & data[[column]] <= value[2]
        },

        # Default
        stop(sprintf("Unknown condition: %s", condition))
      )

      # Handle NAs if requested
      if (na.rm && condition %in% c(">", ">=", "<", "<=", "==", "!=")) {
        result[is.na(result)] <- FALSE
      }

      filter_results[[i]] <- result

    } else {
      stop(sprintf("Invalid filter type at position %d", i))
    }
  }

  # Combine filter results
  if (filter_mode == "all") {
    # AND logic - all filters must pass
    final_filter <- Reduce("&", filter_results)
  } else if (filter_mode == "any") {
    # OR logic - any filter must pass
    final_filter <- Reduce("|", filter_results)
  } else {
    stop("filter_mode must be 'all' or 'any'")
  }

  # Apply final filter
  filtered_data <- data[final_filter, ]

  # Select columns if specified
  if (!is.null(keep_columns)) {
    available_cols <- intersect(keep_columns, names(filtered_data))
    if (length(available_cols) < length(keep_columns)) {
      missing <- setdiff(keep_columns, names(filtered_data))
      warning(sprintf("Columns not found: %s", paste(missing, collapse = ", ")))
    }
    filtered_data <- filtered_data[, available_cols, with = FALSE]
  }

  # Report results if verbose
  if (verbose) {
    cat(sprintf("\nFiltering summary:\n"))
    cat(sprintf("  Original rows: %d\n", original_rows))
    cat(sprintf("  Filtered rows: %d\n", nrow(filtered_data)))
    cat(sprintf("  Removed: %d (%.1f%%)\n",
                original_rows - nrow(filtered_data),
                100 * (original_rows - nrow(filtered_data)) / original_rows))
  }

  return(filtered_data)
}

#' Create filter specifications easily
#'
#' @description
#' Helper function to create filter specifications with validation
#'
#' @param column Column name to filter on
#' @param condition Filter condition (e.g., ">", "==", "%in%", "matches")
#' @param value Value(s) to compare against
#' @param na.rm Remove NAs when filtering (default: FALSE)
#'
#' @return A filter specification list
#'
#' @examples
#' # Create a simple filter
#' f1 <- make_filter("CPM_heart", ">", 100)
#'
#' # Create multiple filters
#' filters <- list(
#'   make_filter("biotype", "==", "protein_coding"),
#'   make_filter("has_eqtl", "==", TRUE),
#'   make_filter("CPM_heart", ">", 50)
#' )
#'
#' @export
make_filter <- function(column, condition, value, na.rm = FALSE) {
  valid_conditions <- c("==", "!=", ">", ">=", "<", "<=",
                       "%in%", "%nin%", "matches", "not_matches",
                       "is_na", "not_na", "between")

  if (!condition %in% valid_conditions) {
    stop(sprintf("Invalid condition '%s'. Must be one of: %s",
                condition, paste(valid_conditions, collapse = ", ")))
  }

  list(
    column = column,
    condition = condition,
    value = value,
    na.rm = na.rm
  )
}

#' Apply common gene filtering presets
#'
#' @description
#' Convenience function with preset filter combinations for common use cases
#'
#' @param data Gene data frame/table
#' @param preset Name of preset filter combination
#' @param custom_filters Additional filters to apply
#' @param ... Additional arguments passed to filter_genes_flexible
#'
#' @return Filtered data
#'
#' @examples
#' # Use a preset
#' filtered <- filter_genes_preset(gene_table, preset = "expressed_coding")
#'
#' # Use preset with custom additions
#' filtered <- filter_genes_preset(
#'   gene_table,
#'   preset = "validated",
#'   custom_filters = list(
#'     make_filter("chr", "==", 3)
#'   )
#' )
#'
#' @export
filter_genes_preset <- function(data,
                               preset = NULL,
                               custom_filters = list(),
                               ...) {

  # Define presets
  presets <- list(
    # Basic filtering for expressed protein-coding genes
    expressed_coding = list(
      make_filter("biotype", "==", "protein_coding"),
      make_filter("CPM_heart", ">", 10, na.rm = TRUE)
    ),

    # Genes ready for validation
    validated = list(
      make_filter("biotype", "==", "protein_coding"),
      make_filter("human_gene_id", "not_na"),
      make_filter("CPM_heart", ">", 50, na.rm = TRUE)
    ),

    # Genes with genetic evidence
    genetic_evidence = list(
      function(df) df$has_eqtl | df$n_variants > 0
    ),

    # High priority candidates
    high_priority = list(
      make_filter("biotype", "==", "protein_coding"),
      make_filter("human_gene_id", "not_na"),
      make_filter("CPM_heart", ">", 100, na.rm = TRUE),
      function(df) df$has_eqtl | df$has_missense
    ),

    # Differentially expressed
    differential = list(
      make_filter("log2FC", ">", 1, na.rm = TRUE),
      make_filter("padj", "<", 0.05, na.rm = TRUE)
    )
  )

  # Get preset filters
  if (!is.null(preset)) {
    if (!preset %in% names(presets)) {
      stop(sprintf("Unknown preset '%s'. Available: %s",
                  preset, paste(names(presets), collapse = ", ")))
    }
    preset_filters <- presets[[preset]]
  } else {
    preset_filters <- list()
  }

  # Combine preset and custom filters
  all_filters <- c(preset_filters, custom_filters)

  # Apply filters
  filter_genes_flexible(data, filters = all_filters, ...)
}

# Helper function for NULL default values (not exported)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}