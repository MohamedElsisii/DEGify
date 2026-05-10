#' Run DEGify bulk RNA-seq downstream analysis
#'
#' DEGify performs automated bulk RNA-seq downstream
#' analysis including differential expression analysis,
#' visualization, heatmap generation, PCA analysis,
#' and Gene Ontology enrichment.
#'
#' @param counts Raw count matrix with genes as rows
#' and samples as columns.
#'
#' @param metadata Sample metadata dataframe containing
#' experimental conditions.
#'
#' @param output_dir Directory for exported results.
#'
#' @param padj_cutoff Adjusted p-value threshold for
#' significant differential expression.
#'
#' @param logfc_cutoff Absolute log2 fold-change threshold
#' for significant differential expression.
#'
#' @param top_n_heatmap Number of top genes displayed
#' in the heatmap.
#'
#' @param top_n_labels Number of top genes labeled
#' in the volcano plot.
#'
#' @param ontology Gene Ontology category for enrichment.
#' One of "BP", "MF", or "CC".
#'
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#' \itemize{
#'   \item DESeq2 dataset object
#'   \item Differential expression results
#'   \item Significant DEGs
#'   \item Volcano plot
#'   \item PCA plot
#'   \item Heatmap matrix
#'   \item GO enrichment results
#' }
#'
#' @examples
#' \dontrun{
#' results <- run_degify(
#'   counts = counts,
#'   metadata = metadata
#' )
#' }
#'
#' @export

run_degify <- function(
    counts,
    metadata,
    output_dir = "DEGify_Output",
    padj_cutoff = 0.05,
    logfc_cutoff = 1,
    top_n_heatmap = 50,
    top_n_labels = 10,
    ontology = "BP",
    seed = 123
) {
  degify_version <- "0.1.0"
  log_msg("Starting DEGify analysis...")
  log_msg(
    paste(
      "Adjusted p-value cutoff:",
      padj_cutoff
    )
  )

  log_msg(
    paste(
      "Log2 fold-change cutoff:",
      logfc_cutoff
    )
  )

  log_msg(
    paste(
      "GO ontology:",
      ontology
    )
  )

  dir.create(
    output_dir,
    showWarnings = FALSE
  )
  output_dir <- normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )

  # -------------------------------------------------------
  # PARAMETER VALIDATION
  # -------------------------------------------------------

  if (
    padj_cutoff <= 0 ||
    padj_cutoff >= 1
  ) {

    stop(
      "padj_cutoff must be between 0 and 1."
    )
  }


  if (logfc_cutoff < 0) {

    stop(
      "logfc_cutoff must be positive."
    )
  }


  if (top_n_heatmap <= 0) {

    stop(
      "top_n_heatmap must be positive."
    )
  }


  if (top_n_labels <= 0) {

    stop(
      "top_n_labels must be positive."
    )
  }

  # -------------------------------------------------------
  # INPUT VALIDATION
  # -------------------------------------------------------

  validate_inputs(
    counts,
    metadata
  )

  # -------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------

  counts <- clean_ensembl_ids(counts)

  counts <- preprocess_counts(counts)


  # -------------------------------------------------------
  # DESeq2
  # -------------------------------------------------------

  deseq_results <- run_deseq(
    counts,
    metadata
  )

  dds <- deseq_results$dds

  res_df <- deseq_results$res_df


  # -------------------------------------------------------
  # Annotation
  # -------------------------------------------------------

  res_df$SYMBOL <- map_ensembl_to_feature(
    res_df$ENSEMBL
  )

  res_df <- map_entrez_ids(res_df)


  # -------------------------------------------------------
  # DEG Processing
  # -------------------------------------------------------

  processed <- process_deg_results(
    res_df = res_df,
    padj_cutoff = padj_cutoff,
    logfc_cutoff = logfc_cutoff
  )

  res_df <- processed$res_df

  sig_deg <- processed$sig_deg


  # -------------------------------------------------------
  # Volcano
  # -------------------------------------------------------

  volcano_plot <- plot_volcano(
    res_df = res_df,
    top_n_labels = top_n_labels,
    output_dir = output_dir
  )


  # -------------------------------------------------------
  # PCA
  # -------------------------------------------------------

  pca_results <- plot_pca(
    dds = dds,
    metadata = metadata,
    output_dir = output_dir
  )

  vsd <- pca_results$vsd

  pca_plot <- pca_results$pca_plot


  # -------------------------------------------------------
  # Heatmap
  # -------------------------------------------------------

  heatmap_matrix <- plot_heatmap(
    vsd,
    res_df,
    metadata,
    top_n_heatmap,
    padj_cutoff,
    output_dir
  )


  # -------------------------------------------------------
  # GO Enrichment
  # -------------------------------------------------------

  go_results <- run_go_enrichment(
    sig_deg = sig_deg,
    ontology = ontology,
    output_dir = output_dir
  )


  # -------------------------------------------------------
  # Export Results
  # -------------------------------------------------------

  export_deg_results(
    res_df,
    sig_deg,
    output_dir
  )

  summary_text <- paste0(

    "DEGify Analysis Summary\n",
    "=======================\n\n",

    "DEGify Version: ",
    degify_version, "\n",

    "GO Ontology: ",
    ontology, "\n\n",

    "Adjusted p-value cutoff: ",
    padj_cutoff, "\n",

    "Log2 fold-change cutoff: ",
    logfc_cutoff, "\n\n",

    "Total genes analyzed: ",
    nrow(res_df), "\n",

    "Significant DEGs: ",
    nrow(sig_deg), "\n",

    "Upregulated genes: ",
    sum(sig_deg$log2FoldChange > 0), "\n",

    "Downregulated genes: ",
    sum(sig_deg$log2FoldChange < 0), "\n"
  )

  writeLines(

    summary_text,

    file.path(
      output_dir,
      "Analysis_Summary.txt"
    )
  )


  writeLines(

    utils::capture.output(utils::sessionInfo()),

    file.path(
      output_dir,
      "sessionInfo.txt"
    )
  )
  generate_report(
    sig_deg = sig_deg,
    go_results = go_results,
    output_dir = output_dir
  )

  log_msg(
    "DEGify analysis completed successfully."
  )


  return(
    list(
      dds = dds,
      res_df = res_df,
      sig_deg = sig_deg,
      volcano_plot = volcano_plot,
      pca_plot = pca_plot,
      heatmap_matrix = heatmap_matrix,
      go_results = go_results
    )
  )
}
