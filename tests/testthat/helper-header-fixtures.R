header_schema <- local({
  repo <- getOption("hbcd.test.root")
  output <- tempfile(fileext = ".tsv")
  status <- system2("uv", c("run", "--project", shQuote(repo), "python",
                           shQuote(file.path(repo, "tests/export_header_columns.py")),
                           shQuote(file.path(repo, "hbcd_headers.json")), shQuote(output)))
  stopifnot(status == 0L)
  schema <- read.delim(output, check.names = FALSE, colClasses = "character",
                       na.strings = character(0))
  unlink(output)
  schema
})

release_columns <- function(kind, variant = 1L) {
  rows <- header_schema[header_schema$kind == kind & header_schema$variant == variant, ]
  rows$column_name[order(as.integer(rows$column_index))]
}

release_frame <- function(kind, rows, variant = 1L) {
  columns <- release_columns(kind, variant)
  stopifnot(all(names(rows) %in% columns))
  result <- as.data.frame(setNames(lapply(columns, function(x) rep(NA, nrow(rows))), columns),
                          check.names = FALSE)
  result[names(rows)] <- rows
  result
}

local_sensitivity_helpers <- function() {
  e <- new.env(parent = globalenv())
  sys.source(file.path(getOption("hbcd.test.root"), "R/sensitivity_covariates.R"), e)
  e
}
