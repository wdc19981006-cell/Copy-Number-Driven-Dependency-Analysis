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
    "V1",
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


