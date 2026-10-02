# ============================================================
# 9. Genome-wide screen
# ============================================================

run_genomewide <- function() {

  folder <- file.path(
    RESULT_ROOT,
    "genomewide_dependency"
  )

  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  cn <- read_gene(
    FILES$cn,
    GENE_A,
    "CN_relative"
  ) %>%
    prepare_cn()

  msg(
    "Reading full CRISPR matrix..."
  )

  dt <- fread(
    FILES$chronos,
    check.names = FALSE,
    showProgress = TRUE
  )

  id_col <- identify_id_column(
    names(dt)
  )

  model_ids <- as.character(
    dt[[id_col]]
  )

  gene_cols <- setdiff(
    names(dt),
    id_col
  )

  genes <- clean_gene(
    gene_cols
  )

  M <- as.matrix(
    dt[
      ,
      ..gene_cols
    ]
  )

  storage.mode(M) <- "double"

  rm(dt)
  gc()

  idx <- match(
    model_ids,
    cn$ModelID
  )

  x <- cn$CN_log[
    idx
  ]

  group <- cn$CN_binary[
    idx
  ]

  valid_model <- is.finite(x)

  M <- M[
    valid_model,
    ,
    drop = FALSE
  ]

  x <- x[
    valid_model
  ]

  group <- group[
    valid_model
  ]

  low_idx <-
    group ==
      "CN-Low"

  nonlow_idx <-
    group ==
      "CN-NonLow"

  n_gene <- ncol(M)

  result <- data.table(

    Gene = genes,

    N = integer(
      n_gene
    ),

    Pearson_r =
      numeric(
        n_gene
      ),

    Pearson_P =
      numeric(
        n_gene
      ),

    N_low =
      integer(
        n_gene
      ),

    N_nonlow =
      integer(
        n_gene
      ),

    Mean_low =
      numeric(
        n_gene
      ),

    Mean_nonlow =
      numeric(
        n_gene
      ),

    Median_low =
      numeric(
        n_gene
      ),

    Median_nonlow =
      numeric(
        n_gene
      ),

    Delta_mean =
      numeric(
        n_gene
      ),

    Delta_median =
      numeric(
        n_gene
      ),

    Wilcoxon_P =
      numeric(
        n_gene
      )

  )

  result[
    ,
    c(
      "Pearson_r",
      "Pearson_P",
      "Mean_low",
      "Mean_nonlow",
      "Median_low",
      "Median_nonlow",
      "Delta_mean",
      "Delta_median",
      "Wilcoxon_P"
    ) :=
      NA_real_
  ]

  for (
    j in seq_len(
      n_gene
    )
  ) {

    y <- M[, j]

    ok <-
      is.finite(x) &
      is.finite(y)

    n <- sum(ok)

    result$N[j] <- n

    if (n >= 10) {

      r <- suppressWarnings(
        cor(
          x[ok],
          y[ok]
        )
      )

      result$Pearson_r[j] <- r

      if (
        is.finite(r) &&
        abs(r) < 1
      ) {

        tval <-
          r *
          sqrt(
            (n - 2) /
              (1 - r^2)
          )

        result$Pearson_P[j] <-
          2 *
          pt(
            -abs(tval),
            df = n - 2
          )

      }

    }

    lv <- y[
      low_idx &
        is.finite(y)
    ]

    nv <- y[
      nonlow_idx &
        is.finite(y)
    ]

    result$N_low[j] <-
      length(lv)

    result$N_nonlow[j] <-
      length(nv)

    if (
      length(lv) >= MIN_N &&
      length(nv) >= MIN_N
    ) {

      result$Mean_low[j] <-
        mean(lv)

      result$Mean_nonlow[j] <-
        mean(nv)

      result$Median_low[j] <-
        median(lv)

      result$Median_nonlow[j] <-
        median(nv)

      result$Delta_mean[j] <-
        mean(lv) -
        mean(nv)

      result$Delta_median[j] <-
        median(lv) -
        median(nv)

      result$Wilcoxon_P[j] <-
        tryCatch(

          wilcox.test(
            lv,
            nv,
            exact = FALSE
          )$p.value,

          error =
            function(e)
              NA_real_

        )

    }

    if (
      j %% 1000 == 0
    ) {

      msg(
        "Completed ",
        j,
        "/",
        n_gene
      )

    }

  }

  result[
    ,
    Pearson_FDR :=
      p.adjust(
        Pearson_P,
        method = "BH"
      )
  ]

  result[
    ,
    Wilcoxon_FDR :=
      p.adjust(
        Wilcoxon_P,
        method = "BH"
      )
  ]

  setorder(
    result,
    Wilcoxon_FDR,
    Delta_median
  )

  result[
    ,
    Rank := .I
  ]

  fwrite(
    result,
    file.path(
      folder,
      "GenomeWide_Dependency.csv"
    )
  )

  hits <- result[
    is.finite(
      Delta_median
    ) &
      is.finite(
        Wilcoxon_FDR
      )
  ]

  hits[
    ,
    negLogFDR :=
      -log10(
        pmax(
          Wilcoxon_FDR,
          1e-300
        )
      )
  ]

  hits[
    ,
    Direction :=
      fifelse(
        Wilcoxon_FDR < 0.05 &
          Delta_median < 0,
        "More dependent in CN-Low",
        fifelse(
          Wilcoxon_FDR < 0.05 &
            Delta_median > 0,
          "More dependent in CN-NonLow",
          "Not significant"
        )
      )
  ]

  top_label <- hits[
    Direction ==
      "More dependent in CN-Low"
  ][
    order(
      Wilcoxon_FDR,
      Delta_median
    )
  ][
    1:min(
      .N,
      8
    ),
    Gene
  ]

  label_genes <- unique(
    c(
      GENE_B,
      top_label
    )
  )

  hits[
    ,
    Label :=
      ifelse(
        Gene %in%
          label_genes,
        Gene,
        NA_character_
      )
  ]

  p <- ggplot(
    hits,
    aes(
      Delta_median,
      negLogFDR
    )
  ) +
    geom_point(
      aes(
        color = Direction
      ),
      alpha = 0.55,
      size = 1.4
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept =
        -log10(0.05),
      linetype = "dashed"
    ) +
    geom_text_repel(
      aes(
        label = Label
      ),
      na.rm = TRUE,
      max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c(
        "More dependent in CN-Low"
        = "#C7473B",
        "More dependent in CN-NonLow"
        = "#3A6EA5",
        "Not significant"
        = "grey70"
      )
    ) +
    labs(
      title =
        paste0(
          GENE_A,
          " CN-Low Genome-wide Dependency Screen"
        ),
      x =
        "Delta median Chronos (CN-Low - CN-NonLow)",
      y =
        "-log10(FDR)",
      color = NULL
    ) +
    theme_classic()

  save_pdf(
    p,
    folder,
    "GenomeWide_Dependency_Volcano",
    width = 8,
    height = 6
  )

  top_hits <- result[
    !is.na(
      Wilcoxon_FDR
    )
  ][
    1:min(
      .N,
      100
    )
  ]

  fwrite(
    top_hits,
    file.path(
      folder,
      "GenomeWide_Dependency_Top100.csv"
    )
  )

  candidate <- result[
    Gene ==
      GENE_B
  ]

  fwrite(
    candidate,
    file.path(
      folder,
      paste0(
        GENE_B,
        "_Candidate_Rank.csv"
      )
    )
  )

  rm(M)
  gc()
}


