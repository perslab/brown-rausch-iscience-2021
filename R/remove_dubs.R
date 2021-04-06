##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export
remove_dubs <- function(srt, pkval) {

  nExp_poi <- 0.075*nrow(srt@meta.data)
  srt <- doubletFinder_v3(srt, PCs = srt@commands$RunUMAP.SCT.pca$dims,
                          sct = TRUE, pN = 0.25, pK = pkval, 
                          nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

}
