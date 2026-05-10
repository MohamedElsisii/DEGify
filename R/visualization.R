#' Generate DEGify volcano plot
#'
#' Creates a volcano plot for differential
#' expression analysis results.
#'
#' @param res_df Differential expression results dataframe.
#'
#' @param top_n_labels Number of top genes labeled
#' on the volcano plot.
#'
#' @param output_dir Directory for exported plots.
#'
#' @return A ggplot2 volcano plot object.
#'
#' @export

# =========================================================
# VOLCANO PLOT
# =========================================================

plot_volcano <- function(
    res_df,
    top_n_labels = 10,
    output_dir = "DEGify_Output"
) {

  log_msg("Generating volcano plot...")


  top_genes <- res_df |>

    subset(
      padj < 0.001 &
        !is.na(SYMBOL) &
        SYMBOL != ""
    ) |>

    dplyr::arrange(padj) |>

    utils::head(top_n_labels)


  volcano_plot <- ggplot2::ggplot(

    res_df,

    ggplot2::aes(
      x = log2FoldChange,
      y = negLog10Padj,
      color = significance
    )

  ) +

    ggplot2::geom_point(
      alpha = 0.7,
      size = 2
    ) +

    ggplot2::scale_color_manual(

      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not Significant" = "gray"
      )
    ) +

    ggrepel::geom_text_repel(

      data = top_genes,

      ggplot2::aes(label = SYMBOL),

      size = 4,

      max.overlaps = 20
    ) +

    ggplot2::theme_minimal() +

    ggplot2::labs(

      title = "Volcano Plot",

      x = "Log2 Fold Change",

      y = "-Log10 Adjusted P-value"
    )


  ggplot2::ggsave(

    file.path(
      output_dir,
      "Volcano_Plot.png"
    ),

    volcano_plot,

    width = 10,

    height = 8,

    dpi = 300
  )

  return(volcano_plot)
}


#' Generate PCA plot
#'
#' Performs variance stabilizing transformation
#' and generates a PCA plot.
#'
#' @param dds DESeq2 dataset object.
#'
#' @param metadata Sample metadata dataframe.
#'
#' @param output_dir Directory for exported plots.
#'
#' @return A list containing:
#' \itemize{
#'   \item Variance stabilized dataset
#'   \item PCA ggplot object
#' }
#'
#' @export

# =========================================================
# PCA PLOT
# =========================================================

plot_pca <- function(
    dds,
    metadata,
    output_dir = "DEGify_Output"
) {

  log_msg("Generating PCA plot...")

  vsd <- DESeq2::varianceStabilizingTransformation(
    dds,
    blind = FALSE
  )

  pca_data <- DESeq2::plotPCA(

    vsd,

    intgroup = "condition",

    returnData = TRUE
  )

  percent_var <- round(
    100 * attr(pca_data, "percentVar")
  )

  pca_plot <- ggplot2::ggplot(

    pca_data,

    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = condition
    )

  ) +

    ggplot2::geom_point(
      size = 5,
      alpha = 0.9
    ) +

    ggplot2::theme_minimal() +

    ggplot2::labs(

      title = "PCA Plot",

      x = paste0(
        "PC1: ",
        percent_var[1],
        "% variance"
      ),

      y = paste0(
        "PC2: ",
        percent_var[2],
        "% variance"
      )
    )


  ggplot2::ggsave(

    file.path(
      output_dir,
      "PCA_Plot.png"
    ),

    pca_plot,

    width = 8,

    height = 6
  )

  return(
    list(
      vsd = vsd,
      pca_plot = pca_plot
    )
  )
}


#' Generate DEGify heatmap
#'
#' Generates a clustered heatmap of top
#' differentially expressed genes.
#'
#' @param vsd Variance stabilized dataset.
#'
#' @param res_df Differential expression results dataframe.
#'
#' @param metadata Sample metadata dataframe.
#'
#' @param top_n_heatmap Number of genes displayed
#' in the heatmap.
#'
#' @param padj_cutoff Adjusted p-value threshold.
#'
#' @param output_dir Directory for exported plots.
#'
#' @return Heatmap expression matrix.
#'
#' @export

# =========================================================
# HEATMAP
# =========================================================

plot_heatmap <- function(
    vsd,
    res_df,
    metadata,
    top_n_heatmap = 50,
    padj_cutoff = 0.05,
    output_dir = "DEGify_Output"
) {

  log_msg("Generating heatmap...")


  top_heatmap_genes <- res_df |>

    dplyr::filter(
      padj < padj_cutoff
    ) |>

    dplyr::arrange(
      padj,
      dplyr::desc(abs(log2FoldChange))
    ) |>

    utils::head(top_n_heatmap)


  heatmap_matrix <- SummarizedExperiment::assay(vsd)[
    top_heatmap_genes$ENSEMBL,
  ]


  gene_labels <- ifelse(
    is.na(top_heatmap_genes$SYMBOL),
    top_heatmap_genes$ENSEMBL,
    top_heatmap_genes$SYMBOL
  )

  rownames(heatmap_matrix) <- make.unique(
    gene_labels
  )


  heatmap_matrix <- t(
    scale(
      t(heatmap_matrix)
    )
  )


  heatmap_matrix <- heatmap_matrix[
    stats::complete.cases(heatmap_matrix),
  ]


  annotation_col <- data.frame(
    condition = metadata$condition
  )

  rownames(annotation_col) <- rownames(metadata)


  condition_levels <- unique(metadata$condition)

  condition_colors <- stats::setNames(
    grDevices::rainbow(length(condition_levels)),
    condition_levels
  )

  annotation_colors <- list(
    condition = condition_colors
  )


  pheatmap::pheatmap(

    heatmap_matrix,

    annotation_col = annotation_col,

    annotation_colors = annotation_colors,

    show_rownames = TRUE,

    show_colnames = FALSE,

    clustering_distance_rows = "euclidean",

    clustering_distance_cols = "euclidean",

    clustering_method = "complete",

    scale = "none",

    fontsize_row = 8,

    main = "Top Differentially Expressed Genes",

    filename = file.path(
      output_dir,
      "Heatmap.png"
    ),

    width = 8,

    height = 10
  )

  return(heatmap_matrix)
}
