
# =========================================================
# RUN DESEQ2
# =========================================================

run_deseq <- function(
    counts,
    metadata,
    design_formula = ~ condition
) {

  log_msg("Running DESeq2 analysis...")

  metadata$condition <- as.factor(
    metadata$condition
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(

    countData = counts,

    colData = metadata,

    design = design_formula
  )

  dds <- DESeq2::DESeq(dds)

  res <- DESeq2::results(dds)

  res_df <- as.data.frame(res)

  res_df$ENSEMBL <- rownames(res_df)

  return(
    list(
      dds = dds,
      res_df = res_df
    )
  )
}


# =========================================================
# PROCESS DEG RESULTS
# =========================================================

process_deg_results <- function(
    res_df,
    padj_cutoff = 0.05,
    logfc_cutoff = 1
) {

  log_msg("Processing DEG results...")

  res_df <- res_df |>

    dplyr::filter(
      !is.na(log2FoldChange),
      !is.na(padj)
    )

  res_df$padj[
    res_df$padj == 0
  ] <- 1e-300


  # -------------------------------------------------------
  # Significance Classification
  # -------------------------------------------------------

  res_df$significance <- dplyr::case_when(

    res_df$padj < padj_cutoff &
      res_df$log2FoldChange > logfc_cutoff
    ~ "Upregulated",

    res_df$padj < padj_cutoff &
      res_df$log2FoldChange < -logfc_cutoff
    ~ "Downregulated",

    TRUE ~ "Not Significant"
  )


  # -------------------------------------------------------
  # Volcano Coordinates
  # -------------------------------------------------------

  res_df$negLog10Padj <- -log10(
    res_df$padj
  )


  # -------------------------------------------------------
  # Significant DEGs
  # -------------------------------------------------------

  sig_deg <- subset(

    res_df,

    padj < padj_cutoff &
      abs(log2FoldChange) > logfc_cutoff
  )

  return(
    list(
      res_df = res_df,
      sig_deg = sig_deg
    )
  )
}
