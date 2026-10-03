# Exercise actual finalizer routing without recalculating a genome-wide screen.
expressions <- parse("scripts/R/adapters.R")
definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name("finalize_outputs")), as.list(expressions))
stopifnot(length(definition)==1L)
e <- new.env(parent=baseenv())
e$PROJECT_ROOT <- "D:/Project"
e$GENE_A <- "GENEB"; e$GENE_B <- "GENEB"
captured <- NULL
e$system2 <- function(command, args) {captured <<- args; 0L}
eval(definition[[1]], e)
for (case in c("GENEA_GENEB", "Custom_Case")) {
  e$RESULT_ROOT <- file.path(e$PROJECT_ROOT,"results",case)
  e$finalize_outputs("genomewide_dependency")
  stopifnot(captured[match("--case",captured)+1L] == shQuote(case))
  stopifnot(captured[match("--candidate",captured)+1L] == "GENEB")
}
cat("PASS: default and custom result folders reach ranking decoration correctly.\n")
