# ============================================================
# Copy Number–Driven Dependency Analysis
# run_analysis.R
#
# R version: 4.5.0
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(broom)
  library(scales)
  library(stringr)
  library(optparse)
})

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

# ============================================================
# 1. CLI
# ============================================================

option_list <- list(

  make_option(
    "--mode",
    type = "character",
    default = "targeted_dependency"
  ),

  make_option(
    "--geneA",
    type = "character",
    default = "ENO1"
  ),

  make_option(
    "--geneB",
    type = "character",
    default = "ENO2"
  ),

  make_option(
    "--project",
    type = "character",
    default = "."
  ),

  make_option(
    "--bootstrap",
    type = "integer",
    default = 1000
  ),

  make_option(
    "--min_group_n",
    type = "integer",
    default = 3
  )

)

opt <- parse_args(
  OptionParser(
    option_list = option_list
  )
)

MODE   <- opt$mode
GENE_A <- toupper(opt$geneA)
GENE_B <- toupper(opt$geneB)

PROJECT_ROOT <- normalizePath(
  opt$project,
  winslash = "/",
  mustWork = TRUE
)

BOOT_R <- opt$bootstrap
MIN_N  <- opt$min_group_n


# ============================================================
# 2. Path
# ============================================================

DEPMAP_DIR <- file.path(
  PROJECT_ROOT,
  "data",
  "raw",
  "depmap",
  "26Q1"
)

RESULT_ROOT <- file.path(
  PROJECT_ROOT,
  "results",
  paste0(
    GENE_A,
    "_",
    GENE_B
  )
)

dir.create(
  RESULT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


FILES <- list(

  model =
    file.path(
      DEPMAP_DIR,
      "Model.csv"
    ),

  condition =
    file.path(
      DEPMAP_DIR,
      "ModelCondition.csv"
    ),

  profiles =
    file.path(
      DEPMAP_DIR,
      "OmicsProfiles.csv"
    ),

  cn =
    file.path(
      DEPMAP_DIR,
      "CopyNumber_WGS_26Q1.csv"
    ),

  expression =
    file.path(
      DEPMAP_DIR,
      "Expression_26Q1.csv"
    ),

  chronos =
    file.path(
      DEPMAP_DIR,
      "CRISPR_Chronos_26Q1.csv"
    ),

  dependency =
    file.path(
      DEPMAP_DIR,
      "CRISPR_GeneDependency_26Q1.csv"
    ),

  damaging =
    file.path(
      DEPMAP_DIR,
      "Mutation_Damaging_26Q1.csv"
    ),

  hotspot =
    file.path(
      DEPMAP_DIR,
      "Mutation_Hotspot_26Q1.csv"
    ),

  signatures =
    file.path(
      DEPMAP_DIR,
      "OmicsSignatures_26Q1.csv"
    ),

  subtype =
    file.path(
      DEPMAP_DIR,
      "MolecularSubtypes_26Q1.csv"
    )

)


# ============================================================
# 3. Helpers
# ============================================================

msg <- function(...) {

  cat(
    "\n",
    paste0(...),
    "\n",
    sep = ""
  )

}


check_file <- function(x) {

  if (!file.exists(x)) {

    stop(
      "Missing file: ",
      x
    )

  }

}


clean_gene <- function(x) {

  x <- sub(
    "\\s*\\([^)]*\\)\\s*$",
    "",
    x
  )

  x <- sub(
    "\\|.*$",
    "",
    x
  )

  x
}


identify_id_column <- function(nms) {

  candidates <- c(
    "ModelID",
    "DepMap_ID",
    "DepMapID",
    "Unnamed: 0",
    "X",
    ""
  )

  for (candidate in candidates) {

    hit <- which(
      nms == candidate
    )

    if (length(hit) > 0) {
      return(nms[hit[1]])
    }

  }

  nms[1]
}


find_gene_col <- function(
    nms,
    gene,
    required = TRUE
) {

  clean <- clean_gene(nms)

  idx <- which(
    toupper(clean) ==
      toupper(gene)
  )

  if (length(idx) > 0) {

    return(
      nms[idx[1]]
    )

  }

  if (required) {

    stop(
      "Gene not found: ",
      gene
    )

  }

  NULL
}


read_gene <- function(
    file,
    gene,
    value_name
) {

  check_file(file)

  header <- fread(
    file,
    nrows = 0,
    check.names = FALSE,
    showProgress = FALSE
  )

  nms <- names(header)

  id_col <- identify_id_column(
    nms
  )

  gene_col <- find_gene_col(
    nms,
    gene
  )

  dat <- fread(
    file,
    select = c(
      id_col,
      gene_col
    ),
    check.names = FALSE,
    showProgress = FALSE
  )

  setnames(
    dat,
    c(
      id_col,
      gene_col
    ),
    c(
      "ModelID",
      value_name
    )
  )

  dat[, ModelID := as.character(ModelID)]

  dat[
    ,
    (value_name) :=
      as.numeric(
        get(value_name)
      )
  ]

  unique(
    dat,
    by = "ModelID"
  )
}


p_star <- function(p) {

  if (is.na(p)) return("NA")
  if (p < 1e-4) return("****")
  if (p < 1e-3) return("***")
  if (p < 1e-2) return("**")
  if (p < 0.05) return("*")

  "ns"
}


p_text <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 2e-16) {
    return("<2e-16")
  }

  format.pval(
    p,
    digits = 3
  )
}


save_pdf <- function(
    plot,
    folder,
    name,
    width = 7,
    height = 6
) {

  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )

  ggsave(
    filename =
      file.path(
        folder,
        paste0(
          name,
          ".pdf"
        )
      ),
    plot = plot,
    width = width,
    height = height
  )
}


# ============================================================
# 4. CN classification
# ============================================================

prepare_cn <- function(dat) {

  dat %>%
    mutate(

      CN_relative =
        as.numeric(
          CN_relative
        ),

      CN_log =
        log2(
          CN_relative + 1
        ),

      CN_binary =
        if_else(
          CN_log < 0.585,
          "CN-Low",
          "CN-NonLow"
        ),

      CN_status =
        case_when(

          CN_log < 0.35
          ~ "Deep CN Loss",

          CN_log < 0.585
          ~ "Shallow CN Loss",

          TRUE
          ~ "CN Non-Low"

        )

    )
}


# ============================================================
# 5. Metadata
# ============================================================

load_model <- function() {

  check_file(
    FILES$model
  )

  x <- fread(
    FILES$model
  )

  if (!"ModelID" %in% names(x)) {
    stop("Model.csv missing ModelID")
  }

  x
}


# ============================================================
# 6. Common pair
# ============================================================

load_target_pair <- function(
    geneA,
    geneB
) {

  cn <- read_gene(
    FILES$cn,
    geneA,
    "CN_relative"
  )

  cn <- prepare_cn(
    cn
  )

  dep <- read_gene(
    FILES$chronos,
    geneB,
    "Chronos"
  )

  inner_join(
    cn,
    dep,
    by = "ModelID"
  ) %>%
    filter(
      is.finite(CN_relative),
      is.finite(CN_log),
      is.finite(Chronos)
    )
}


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


# ============================================================
# 12. Reverse
# ============================================================

run_reverse <- function() {

  old_A <- GENE_A
  old_B <- GENE_B

  assign(
    "GENE_A",
    old_B,
    envir = .GlobalEnv
  )

  assign(
    "GENE_B",
    old_A,
    envir = .GlobalEnv
  )

  assign(
    "RESULT_ROOT",
    file.path(
      PROJECT_ROOT,
      "results",
      paste0(
        old_B,
        "_",
        old_A
      )
    ),
    envir = .GlobalEnv
  )

  dir.create(
    RESULT_ROOT,
    recursive = TRUE,
    showWarnings = FALSE
  )

  run_targeted()

  assign(
    "GENE_A",
    old_A,
    envir = .GlobalEnv
  )

  assign(
    "GENE_B",
    old_B,
    envir = .GlobalEnv
  )
}


# ============================================================
# 13. Dispatcher
# ============================================================

msg(
  "R version: ",
  R.version.string
)

msg(
  "Mode: ",
  MODE
)

msg(
  "Gene A: ",
  GENE_A
)

msg(
  "Gene B: ",
  GENE_B
)


if (
  MODE ==
    "depmap_cn_expression"
) {

  run_cn_expression()

} else if (
  MODE ==
    "targeted_dependency"
) {

  run_targeted()

} else if (
  MODE ==
    "genomewide_dependency"
) {

  run_genomewide()

} else if (
  MODE ==
    "lineage_dependency"
) {

  run_lineage()

} else if (
  MODE ==
    "adjusted_dependency"
) {

  run_adjusted()

} else if (
  MODE ==
    "reverse_dependency"
) {

  run_reverse()

} else if (
  MODE ==
    "full"
) {

  run_cn_expression()

  gc()

  run_genomewide()

  gc()

  run_targeted()

  gc()

  run_lineage()

  gc()

  run_adjusted()

  gc()

  run_reverse()

  gc()

} else {

  stop(
    "Unsupported mode: ",
    MODE
  )

}


msg(
  "Analysis completed."
)