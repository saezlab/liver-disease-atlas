ma_qc = function(eset_raw, nuse_spread = 0.1, rle_spread = 0.1) {
  
  plmfit = fitProbeLevelModel(eset_raw)
  
  discard = F
  
  nuse  = apply(oligo::NUSE(plmfit, type = "values"), 2, median, na.rm = T)
  discard = discard | nuse > 1 + nuse_spread | nuse < 1 - nuse_spread
  
  rle = apply(oligo::RLE(plmfit, type="value"), 2, median, na.rm = T)
  discard = discard | rle > rle_spread | rle < -rle_spread
  
  if (any(discard)) {
    warning(paste("Discarding in total", 
                  sum(discard), 
                  "arrays:", 
                  paste(names(discard[discard == T]), collapse = " ")))
    eset_raw = eset_raw[,!discard]
  } else {
    message("All arrays passed quality control.")
  }
  
  return(eset_raw)
}


ma_annotate = function(eset, platforms) {
  mat = exprs(eset)
  annotation = platforms[[eset@annotation]]
  rownames(mat) = annotate::getSYMBOL(as.vector(rownames(mat)), annotation)
  expr = limma::avereps(mat[!is.na(rownames(mat)),])
  
  return(expr)
}