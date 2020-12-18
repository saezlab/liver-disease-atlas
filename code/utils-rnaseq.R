preprocess_count_matrix <- function(count_matrix, verbose = TRUE) {

  # removes genes with constant expression across all samples
  # (including genes with 0 counts)
  keep <- apply(count_matrix, 1, var) != 0

  res <- log2(count_matrix[keep, ] + 1)

  if (verbose) {
    message("Discarding ", sum(!keep), " genes \nKeeping ", sum(keep), " genes")
  }
  return(res)
}

voom_normalization <- function(count_matrix) {

  # filter low read counts, TMM normalization and logCPM transformation
  keep <- filterByExpr(count_matrix)

  norm <- count_matrix[keep, , keep.lib.sizes = F] %>%
    calcNormFactors() %>%
    voom() %>%
    pluck("E")

  message("Discarding ", sum(!keep), " genes \nKeeping ", sum(keep), " genes")

  return(norm)
}
