##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param psrt processed seurat object

##' @return
##' @author dylanmr
##' @export

calc_pK <- function(psrt) {

  sweep.res.list <- DoubletFinder::paramSweep_v3(seu = psrt, 
                                                 PCs = psrt@commands$RunUMAP.SCT.pca$dims, 
                                                 sct=T)
  
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

}
