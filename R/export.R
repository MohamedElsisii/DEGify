# =========================================================
# EXPORT RESULTS
# =========================================================

export_deg_results <- function(
    res_df,
    sig_deg,
    output_dir = "DEGify_Output"
) {

  log_msg("Exporting results...")


  res_df <- res_df |>

    dplyr::select(
      ENSEMBL,
      SYMBOL,
      ENTREZID,
      significance,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      negLog10Padj
    )


  sig_deg <- sig_deg |>

    dplyr::select(
      ENSEMBL,
      SYMBOL,
      ENTREZID,
      significance,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      negLog10Padj
    ) |>

    dplyr::arrange(padj)


  utils::write.csv(

    res_df,

    file.path(
      output_dir,
      "All_DEGs.csv"
    )
  )


  utils::write.csv(

    sig_deg,

    file.path(
      output_dir,
      "Significant_DEGs.csv"
    )
  )
}
