# =========================================================
# MAP ENSEMBL TO SYMBOL
# =========================================================

map_ensembl_to_feature <- function(
    ensembl_vec,
    mart = NULL
) {

  ens_clean <- sub(
    "\\..*$",
    "",
    ensembl_vec
  )

  mapping <- NULL


  # -------------------------------------------------------
  # biomaRt Mapping
  # -------------------------------------------------------

  if (!is.null(mart)) {

    log_msg(
      "Mapping Ensembl -> SYMBOL via biomaRt..."
    )

    mapping <- tryCatch({

      biomaRt::getBM(

        attributes = c(
          "ensembl_gene_id",
          "external_gene_name"
        ),

        filters = "ensembl_gene_id",

        values = unique(ens_clean),

        mart = mart

      ) |>

        tibble::as_tibble() |>

        dplyr::rename(
          ENSEMBL = ensembl_gene_id,
          SYMBOL = external_gene_name
        ) |>

        dplyr::mutate(
          SYMBOL = ifelse(
            SYMBOL == "",
            NA,
            SYMBOL
          )
        )

    }, error = function(e) {

      log_msg(
        "biomaRt mapping failed. Falling back to org.Hs.eg.db"
      )

      NULL
    })
  }


  # -------------------------------------------------------
  # Offline Mapping
  # -------------------------------------------------------

  if (is.null(mapping)) {

    log_msg(
      "Mapping Ensembl -> SYMBOL using org.Hs.eg.db..."
    )

    mapping <- AnnotationDbi::select(

      x = org.Hs.eg.db::org.Hs.eg.db,

      keys = unique(ens_clean),

      keytype = "ENSEMBL",

      columns = c(
        "SYMBOL",
        "ENTREZID"
      )

    ) |>

      tibble::as_tibble() |>

      dplyr::mutate(
        SYMBOL = ifelse(
          SYMBOL == "",
          NA,
          SYMBOL
        )
      ) |>

      dplyr::mutate(
        FEATURE = ifelse(
          !is.na(SYMBOL),
          SYMBOL,
          ENSEMBL
        )
      ) |>

      dplyr::select(
        ENSEMBL,
        FEATURE
      )

  } else {

    mapping <- mapping |>

      dplyr::mutate(
        FEATURE = SYMBOL
      ) |>

      dplyr::select(
        ENSEMBL,
        FEATURE
      )
  }


  # -------------------------------------------------------
  # Preserve Order
  # -------------------------------------------------------

  id_map <- mapping |>
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)

  feat <- id_map$FEATURE[
    match(
      ens_clean,
      id_map$ENSEMBL
    )
  ]

  return(feat)
}


# =========================================================
# MAP ENTREZ IDS
# =========================================================

map_entrez_ids <- function(res_df) {

  log_msg("Mapping ENTREZ IDs...")

  entrez_map <- AnnotationDbi::select(

    org.Hs.eg.db::org.Hs.eg.db,

    keys = res_df$ENSEMBL,

    columns = "ENTREZID",

    keytype = "ENSEMBL"
  )

  entrez_map <- entrez_map |>
    dplyr::distinct(ENSEMBL, .keep_all = TRUE)

  res_df$ENTREZID <- entrez_map$ENTREZID[
    match(
      res_df$ENSEMBL,
      entrez_map$ENSEMBL
    )
  ]

  return(res_df)
}
