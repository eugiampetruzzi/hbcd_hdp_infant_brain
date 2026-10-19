# Run with: uv run Rscript tests/testthat.R [--require-combat]
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% "--require-combat")) {
  stop("Usage: Rscript tests/testthat.R [--require-combat]", call. = FALSE)
}
required <- c("testthat", "withr", "mgcv")
if ("--require-combat" %in% args) {
  required <- c(required, "ComBatFamily")
}
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing R test dependencies: ", paste(missing, collapse = ", "), call. = FALSE)
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- normalizePath(file.path(dirname(script), ".."))
options(hbcd.test.root = root)
testthat::test_dir(file.path(root, "tests/testthat"), reporter = "summary",
                   stop_on_failure = TRUE)
