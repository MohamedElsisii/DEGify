# =========================================================
# CLEAN ENSEMBL IDS
# =========================================================

clean_ensembl_ids <- function(counts) {

  log_msg("Cleaning ENSEMBL IDs...")

  clean_ids <- sub(
    "\\..*",
    "",
    rownames(counts)
  )

  counts <- as.data.frame(counts)

  counts$ENSEMBL <- clean_ids

  counts <- counts |>
    dplyr::group_by(ENSEMBL) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), sum))

  counts <- as.data.frame(counts)

  rownames(counts) <- counts$ENSEMBL

  counts$ENSEMBL <- NULL

  return(counts)
}


# =========================================================
# PREPROCESS COUNTS
# =========================================================

preprocess_counts <- function(
    counts,
    min_count = 10,
    min_samples = 2
) {

  log_msg("Preprocessing counts...")

  counts <- round(as.matrix(counts))

  mode(counts) <- "integer"

  keep <- rowSums(counts >= min_count) >= min_samples

  counts <- counts[keep, ]

  return(counts)
}

# -------------------------------------------------------
# Validation
# -------------------------------------------------------

validate_inputs <- function(
    counts,
    metadata
) {

  log_msg("Validating inputs...")


  # -----------------------------------------------------
  # Count Matrix Validation
  # -----------------------------------------------------

  if (
    !is.matrix(counts) &&
    !is.data.frame(counts)
  ) {

    stop(
      "counts must be a matrix or data.frame."
    )
  }


  # -----------------------------------------------------
  # Metadata Validation
  # -----------------------------------------------------

  if (!is.data.frame(metadata)) {

    stop(
      "metadata must be a data.frame."
    )
  }


  # -----------------------------------------------------
  # Condition Column Validation
  # -----------------------------------------------------

  if (
    !"condition" %in% colnames(metadata)
  ) {

    stop(
      "metadata must contain a 'condition' column."
    )
  }


  # -----------------------------------------------------
  # Sample Number Validation
  # -----------------------------------------------------

  if (ncol(counts) < 2) {

    stop(
      "At least two samples are required."
    )
  }


  # -----------------------------------------------------
  # Sample Matching Validation
  # -----------------------------------------------------

  if (
    !all(
      colnames(counts) ==
      rownames(metadata)
    )
  ) {

    stop(
      "Sample names between counts and metadata do not match."
    )
  }


  return(TRUE)
}
