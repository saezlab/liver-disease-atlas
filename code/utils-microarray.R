# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

#' This function performs a quality check for arrays based on the relative log
#' expression values (RLE) and the normalized unscaled standard errors (NUSE).
#'
#' @param eset_raw GeneFeatureSet constructed by oligo::read.celfiles().
#' @param nuse_spread numeric; allowed NUSE spread centered around 1.
#' @param rle_spread numeric; allowed RLE spread centered around 0.
#'
#' @return Provided GeneFeatureSet with discarded samples that did not pass the
#' quality check.
ma_qc <- function(eset_raw, nuse_spread = 0.1, rle_spread = 0.1) {
  plmfit <- fitProbeLevelModel(eset_raw)

  discard <- F

  nuse <- apply(oligo::NUSE(plmfit, type = "values"), 2, median, na.rm = T)
  discard <- discard | nuse > 1 + nuse_spread | nuse < 1 - nuse_spread

  rle <- apply(oligo::RLE(plmfit, type = "value"), 2, median, na.rm = T)
  discard <- discard | rle > rle_spread | rle < -rle_spread

  if (any(discard)) {
    warning(paste(
      "Discarding in total",
      sum(discard),
      "arrays:",
      paste(names(discard[discard == T]), collapse = " ")
    ))
    eset_raw <- eset_raw[, !discard]
  } else {
    message("All arrays passed quality control.")
  }

  return(eset_raw)
}

#' Function to annotate and summarize probe ID of various microarray platforms.
#'
#' @param eset ExpressionSet with probe IDs as features/rownames.
#' @param platforms named list matching platform design to annotation package.
#' e.g. list(pd.hg.u133a = "hgu133a.db")
#'
#' @return ExpressionSet with annotated gene symbols. Duplicated gene symbols
#' were averaged.
ma_annotate <- function(eset, platforms) {
  mat <- exprs(eset)
  annotation <- platforms[[eset@annotation]]
  rownames(mat) <- annotate::getSYMBOL(as.vector(rownames(mat)), annotation)
  expr <- limma::avereps(mat[!is.na(rownames(mat)), ])

  return(expr)
}
