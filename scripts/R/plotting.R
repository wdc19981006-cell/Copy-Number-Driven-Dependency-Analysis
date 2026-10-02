# ============================================================
# 7. CN -> Expression
# ============================================================

run_cn_expression <- function() {

  folder <- file.path(
    RESULT_ROOT,
    "cn_expression"
  )

  cn <- read_gene(
    FILES$cn,
    GENE_A,
    "CN_relative"
  ) %>%
    prepare_cn()

  exprA <- read_gene(
    FILES$expression,
    GENE_A,
    "Expression_A"
  )

  depB <- read_gene(
    FILES$chronos,
    GENE_B,
    "Chronos_B"
  )

  a <- inner_join(
    cn,
    exprA,
    by = "ModelID"
  ) %>%
    filter(
      is.finite(CN_relative),
      is.finite(Expression_A)
    )

  b <- inner_join(
    exprA,
    depB,
    by = "ModelID"
  ) %>%
    filter(
      is.finite(Expression_A),
      is.finite(Chronos_B)
    )

  p1_test <- cor.test(
    a$CN_relative,
    a$Expression_A,
    method = "pearson"
  )

  s1_test <- cor.test(
    a$CN_relative,
    a$Expression_A,
    method = "spearman",
    exact = FALSE
  )

  p2_test <- cor.test(
    b$Expression_A,
    b$Chronos_B,
    method = "pearson"
  )

  s2_test <- cor.test(
    b$Expression_A,
    b$Chronos_B,
    method = "spearman",
    exact = FALSE
  )

  p1 <- ggplot(
    a,
    aes(
      CN_relative,
      Expression_A
    )
  ) +
    geom_point(
      alpha = 0.35,
      size = 1.2
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black"
    ) +
    labs(
      title =
        paste0(
          GENE_A,
          " Copy Number vs Expression"
        ),
      x =
        paste0(
          GENE_A,
          " Relative Copy Number"
        ),
      y =
        paste0(
          GENE_A,
          " Expression"
        )
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.15,
      label = paste0(
        "Pearson r = ",
        round(
          p1_test$estimate,
          3
        ),
        "\nSpearman rho = ",
        round(
          s1_test$estimate,
          3
        ),
        "\nN = ",
        nrow(a)
      )
    ) +
    theme_classic()

  p2 <- ggplot(
    b,
    aes(
      Expression_A,
      Chronos_B
    )
  ) +
    geom_point(
      alpha = 0.35,
      size = 1.2
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black"
    ) +
    labs(
      title =
        paste0(
          GENE_A,
          " Expression vs ",
          GENE_B,
          " Dependency"
        ),
      x =
        paste0(
          GENE_A,
          " Expression"
        ),
      y =
        paste0(
          GENE_B,
          " Chronos"
        )
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.15,
      label = paste0(
        "Pearson r = ",
        round(
          p2_test$estimate,
          3
        ),
        "\nSpearman rho = ",
        round(
          s2_test$estimate,
          3
        ),
        "\nN = ",
        nrow(b)
      )
    ) +
    theme_classic()

  save_pdf(
    p1 + p2,
    folder,
    "CN_Expression_Dependency",
    width = 12,
    height = 5
  )

  stats <- data.frame(

    analysis = c(
      "CN_vs_expression",
      "expression_vs_dependency"
    ),

    pearson = c(
      p1_test$estimate,
      p2_test$estimate
    ),

    pearson_p = c(
      p1_test$p.value,
      p2_test$p.value
    ),

    spearman = c(
      s1_test$estimate,
      s2_test$estimate
    ),

    spearman_p = c(
      s1_test$p.value,
      s2_test$p.value
    )

  )

  fwrite(
    stats,
    file.path(
      folder,
      "CN_Expression_Statistics.csv"
    )
  )
}



# ============================================================
# 8. Targeted dependency
# ============================================================

run_targeted <- function() {

  folder <- file.path(
    RESULT_ROOT,
    "targeted_dependency"
  )

  pair <- load_target_pair(
    GENE_A,
    GENE_B
  )

  pearson <- cor.test(
    pair$CN_log,
    pair$Chronos,
    method = "pearson"
  )

  spearman <- cor.test(
    pair$CN_log,
    pair$Chronos,
    method = "spearman",
    exact = FALSE
  )

  pair$CN_binary <- factor(
    pair$CN_binary,
    levels = c(
      "CN-NonLow",
      "CN-Low"
    )
  )

  wt <- wilcox.test(
    Chronos ~ CN_binary,
    data = pair,
    exact = FALSE
  )

  yr <- diff(
    range(
      pair$Chronos,
      na.rm = TRUE
    )
  )

  y_pos <-
    max(
      pair$Chronos,
      na.rm = TRUE
    ) +
    yr * 0.12

  bracket <- data.frame(
    group1 = "CN-NonLow",
    group2 = "CN-Low",
    y.position = y_pos,
    label = paste0(
      p_star(
        wt$p.value
      ),
      "\nP ",
      p_text(
        wt$p.value
      )
    )
  )

  p_scatter <- ggplot(
    pair,
    aes(
      CN_log,
      Chronos
    )
  ) +
    geom_point(
      alpha = 0.4
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "black"
    ) +
    labs(
      title =
        paste0(
          GENE_A,
          " CN vs ",
          GENE_B,
          " Dependency"
        ),
      x =
        paste0(
          GENE_A,
          " log2(relative CN + 1)"
        ),
      y =
        paste0(
          GENE_B,
          " Chronos"
        )
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.05,
      vjust = 1.15,
      label = paste0(
        "Pearson r = ",
        round(
          pearson$estimate,
          3
        ),
        "\nSpearman rho = ",
        round(
          spearman$estimate,
          3
        ),
        "\nN = ",
        nrow(pair)
      )
    ) +
    theme_classic()

  p_box <- ggplot(
    pair,
    aes(
      CN_binary,
      Chronos,
      fill = CN_binary
    )
  ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.25,
      size = 1
    ) +
    stat_pvalue_manual(
      bracket,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      bracket.size = 0.7,
      tip.length = 0.02
    ) +
    expand_limits(
      y = y_pos + yr * 0.10
    ) +
    scale_fill_manual(
      values = c(
        "CN-NonLow" = "grey65",
        "CN-Low" = "#C7473B"
      )
    ) +
    labs(
      title =
        paste0(
          GENE_B,
          " Dependency by ",
          GENE_A,
          " CN Status"
        ),
      x = NULL,
      y =
        paste0(
          GENE_B,
          " Chronos"
        )
    ) +
    theme_classic() +
    theme(
      legend.position = "none"
    )

  save_pdf(
    p_scatter + p_box,
    folder,
    "Targeted_Dependency",
    width = 11,
    height = 5.5
  )

  summary <- pair %>%
    group_by(
      CN_binary
    ) %>%
    summarise(
      N = n(),
      mean_chronos =
        mean(
          Chronos
        ),
      median_chronos =
        median(
          Chronos
        ),
      sd =
        sd(
          Chronos
        ),
      .groups = "drop"
    )

  fwrite(
    summary,
    file.path(
      folder,
      "Group_Summary.csv"
    )
  )

  stats <- data.frame(
    geneA = GENE_A,
    geneB = GENE_B,
    N = nrow(pair),
    pearson_r =
      unname(
        pearson$estimate
      ),
    pearson_p =
      pearson$p.value,
    spearman_rho =
      unname(
        spearman$estimate
      ),
    spearman_p =
      spearman$p.value,
    wilcoxon_p =
      wt$p.value
  )

  fwrite(
    stats,
    file.path(
      folder,
      "Statistics.csv"
    )
  )

  fwrite(
    pair,
    file.path(
      folder,
      "CellLines.csv"
    )
  )
}


