generate_report <- function(
    sig_deg,
    go_results,
    output_dir = "DEGify_Output"
) {

  log_msg("Generating HTML report...")

  if (!dir.exists(output_dir)) {

    dir.create(
      output_dir,
      recursive = TRUE
    )
  }

  report_template <- system.file(
    "report_template.Rmd",
    package = "DEGify"
  )

  stopifnot(
    file.exists(file.path(output_dir, "Volcano_Plot.png")),
    file.exists(file.path(output_dir, "PCA_Plot.png")),
    file.exists(file.path(output_dir, "Heatmap.png")),
    file.exists(file.path(output_dir, "GO_Barplot.png"))
  )
  Sys.sleep(2)

  rmarkdown::render(

    input = report_template,

    output_file = "DEGify_Report.html",

    output_dir = output_dir,

    params = list(

      sig_deg = sig_deg,

      go_results = as.data.frame(
        go_results$go_results
      ),

      volcano_path = file.path(
        output_dir,
        "Volcano_Plot.png"
      ),

      pca_path = file.path(
        output_dir,
        "PCA_Plot.png"
      ),

      heatmap_path = file.path(
        output_dir,
        "Heatmap.png"
      ),

      go_barplot_path = file.path(
        output_dir,
        "GO_Barplot.png"
      )
    ),

    envir = new.env(parent = globalenv())
  )
}
