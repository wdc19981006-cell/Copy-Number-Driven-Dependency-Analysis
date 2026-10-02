# ============================================================
# 10. Lineage forest
# ============================================================

run_lineage <- function() {

  folder <- file.path(
    RESULT_ROOT,
    "lineage_dependency"
  )

  pair <- load_target_pair(
    GENE_A,
    GENE_B
  )

  model <- load_model()

  dat <- inner_join(
    pair,
    model,
    by = "ModelID"
  ) %>%
    filter(
      !is.na(
        OncotreeLineage
      )
    )

  lineages <- split(
    dat,
    dat$OncotreeLineage
  )

  out <- list()

  set.seed(1234)

  k <- 1

  for (
    nm in names(
      lineages
    )
  ) {

    d <- lineages[[nm]]

    low <- d$Chronos[
      d$CN_binary ==
        "CN-Low"
    ]

    nonlow <- d$Chronos[
      d$CN_binary ==
        "CN-NonLow"
    ]

    low <- low[
      is.finite(low)
    ]

    nonlow <- nonlow[
      is.finite(nonlow)
    ]

    if (
      length(low) < MIN_N ||
      length(nonlow) < MIN_N
    ) {
      next
    }

    delta <-
      median(low) -
      median(nonlow)

    boot <- replicate(
      BOOT_R,
      median(
        sample(
          low,
          replace = TRUE
        )
      ) -
        median(
          sample(
            nonlow,
            replace = TRUE
          )
        )
    )

    ci <- quantile(
      boot,
      c(
        0.025,
        0.975
      ),
      na.rm = TRUE
    )

    p <- wilcox.test(
      low,
      nonlow,
      exact = FALSE
    )$p.value

    out[[k]] <- data.frame(

      Lineage = nm,

      N_low =
        length(low),

      N_nonlow =
        length(nonlow),

      Median_low =
        median(low),

      Median_nonlow =
        median(nonlow),

      Delta_median =
        delta,

      CI_low =
        ci[1],

      CI_high =
        ci[2],

      P =
        p

    )

    k <- k + 1

  }

  if (
    length(out) == 0
  ) {

    warning(
      "No lineage had sufficient samples."
    )

    return(
      invisible(NULL)
    )
  }

  res <- bind_rows(
    out
  )

  res$FDR <- p.adjust(
    res$P,
    method = "BH"
  )

  res <- res %>%
    arrange(
      Delta_median
    )

  fwrite(
    res,
    file.path(
      folder,
      "Lineage_Dependency.csv"
    )
  )

  res$Lineage <- factor(
    res$Lineage,
    levels = rev(
      res$Lineage
    )
  )

  p <- ggplot(
    res,
    aes(
      Delta_median,
      Lineage
    )
  ) +
    geom_segment(
      aes(
        x = CI_low,
        xend = CI_high,
        yend = Lineage
      )
    ) +
    geom_point(
      aes(
        shape =
          FDR < 0.05
      ),
      size = 3
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    scale_shape_manual(
      values =
        c(
          "TRUE" = 16,
          "FALSE" = 1
        )
    ) +
    labs(
      title =
        paste0(
          GENE_B,
          " Dependency in ",
          GENE_A,
          " CN-Low Models"
        ),
      x =
        "Delta median Chronos (CN-Low - CN-NonLow)",
      y = NULL,
      shape = "FDR < 0.05"
    ) +
    theme_classic()

  save_pdf(
    p,
    folder,
    "Lineage_Dependency_Forest",
    width = 8,
    height =
      max(
        5,
        0.32 *
          nrow(res) +
          2
      )
  )
}



# ============================================================
# 11. Adjusted regression
# ============================================================

run_adjusted <- function() {

  folder <- file.path(
    RESULT_ROOT,
    "adjusted_dependency"
  )

  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  pair <- load_target_pair(
    GENE_A,
    GENE_B
  )

  model <- load_model()

  dat <- inner_join(
    pair,
    model,
    by = "ModelID"
  ) %>%
    filter(
      !is.na(
        OncotreeLineage
      )
    )

  dat$OncotreeLineage <-
    factor(
      dat$OncotreeLineage
    )

  fit1 <- lm(
    Chronos ~
      CN_log +
      OncotreeLineage,
    data = dat
  )

  fit2 <- lm(
    Chronos ~
      I(
        CN_binary ==
          "CN-Low"
      ) +
      OncotreeLineage,
    data = dat
  )

  fwrite(
    broom::tidy(
      fit1
    ),
    file.path(
      folder,
      "Continuous_CN_Adjusted.csv"
    )
  )

  fwrite(
    broom::tidy(
      fit2
    ),
    file.path(
      folder,
      "CNLow_Adjusted.csv"
    )
  )
}


