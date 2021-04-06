##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param glia_cluster
##' @return
##' @author dylanmr
##' @export
runsil <- function(merged_seur, ...) {

  #updating to run
  pc <- merged_seur@reductions$pca@cell.embeddings[,1:30]
  distance <- parallelDist::parDist(pc, method="euclidean")
  sil <- compute.sil(x = Idents(merged_seur), dist = distance)
  
  merged_seur$silhouette <- sil[,1]
  merged_seur$remove <- "No"
  merged_seur$remove[which(merged_seur$silhouette < 0)] <- "Yes"
  
  # Remove cells with low silhouette coefficient and recluster
  merged_seur <- subset(merged_seur, subset = remove == "Yes", invert = T)
  merged_seur <- quickclus(merged_seur)
  Idents(merged_seur) <- "orig_labels"
  return(merged_seur)
}
