# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

#' Function to preprocess a RNA-seq count matrix.
#'
#' @description Genes with constant expression are removed and counts are
#' transformed to log2 space.
#'
#' @param count_matrix count matrix with genes in rows and samples in columns.
#' @param verbose logical indiciating whether it should be printed how many
#' genes are discarded.
#'
#' @return matrix with only non-constant genes containing counts in log2 space.
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


#' Function to perform RNA-seq normalization
#'
#' @param count_matrix count matrix with genes in rows and samples in columns.
#'
#' @return Matrix with normalized expression data.
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
