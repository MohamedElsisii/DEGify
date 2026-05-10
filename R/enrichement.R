#' Run Gene Ontology enrichment analysis
#'
#' Performs GO enrichment analysis on
#' significant differentially expressed genes.
#'
#' @param sig_deg Significant DEG dataframe.
#'
#' @param ontology GO ontology category.
#' One of "BP", "MF", or "CC".
#'
#' @param output_dir Directory for exported results.
#'
#' @return A list containing:
#' \itemize{
#'   \item GO enrichment results
#'   \item GO barplot
#'   \item GO dotplot
#' }
#'
#' @export

# =========================================================
# GO ENRICHMENT
# =========================================================

run_go_enrichment <- function(
    sig_deg,
    ontology = "BP",
    output_dir = "DEGify_Output"
) {

  # -------------------------------------------------------
  # ONTOLOGY VALIDATION
  # -------------------------------------------------------

  valid_ontologies <- c(
    "BP",
    "MF",
    "CC"
  )

  if (!ontology %in% valid_ontologies) {

    stop(
      "ontology must be one of: BP, MF, CC"
    )
  }


  log_msg("Running GO enrichment...")


  sig_entrez <- sig_deg$ENTREZID

  sig_entrez <- sig_entrez[
    !is.na(sig_entrez)
  ]


  if (length(sig_entrez) == 0) {

    stop(
      "No ENTREZ IDs available for enrichment."
    )
  }


  go_results <- clusterProfiler::enrichGO(

    gene = sig_entrez,

    OrgDb = org.Hs.eg.db::org.Hs.eg.db,

    keyType = "ENTREZID",

    ont = ontology,

    pAdjustMethod = "BH",

    pvalueCutoff = 0.05,

    qvalueCutoff = 0.05,

    readable = TRUE
  )


  if (
    is.null(go_results) ||
    nrow(as.data.frame(go_results)) == 0
  ) {

    stop(
      "No enriched GO terms found."
    )
  }


  utils::write.csv(

    as.data.frame(go_results),

    file.path(
      output_dir,
      "GO_Enrichment.csv"
    )
  )


  go_barplot <- graphics::barplot(

    go_results,

    showCategory = 10,

    title = paste(
      "GO",
      ontology,
      "Enrichment"
    )
  )


  ggplot2::ggsave(

    file.path(
      output_dir,
      "GO_Barplot.png"
    ),

    go_barplot,

    width = 10,

    height = 7
  )


  go_dotplot <- enrichplot::dotplot(

    go_results,

    showCategory = 10,

    title = paste(
      "GO",
      ontology,
      "Enrichment"
    )
  )


  ggplot2::ggsave(

    file.path(
      output_dir,
      "GO_Dotplot.png"
    ),

    go_dotplot,

    width = 10,

    height = 7
  )


  return(
    list(
      go_results = go_results,
      go_barplot = go_barplot,
      go_dotplot = go_dotplot
    )
  )
}
