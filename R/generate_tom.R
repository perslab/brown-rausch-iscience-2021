##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param sft soft power calculation 
##' @param dat raw exp data
##' @return
##' @author dylanmr
##' @export

generate_tom <- function(sft, dat, networkType) {
  
  adj <- adjacency(dat, type = "signed", power = sft)
  diag(adj) <- 0
  TOM <- TOMsimilarityFromExpr(dat, networkType = networkType, TOMType = "signed", power = sft, maxPOutliers = 0.05)
  colnames(TOM) <- rownames(TOM) <- colnames(dat)
  dissTOM <- 1 - TOM
  return(dissTOM)
  
}
