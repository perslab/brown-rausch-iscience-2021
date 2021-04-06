##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seurobj seurat object
##' @param clusid which metadata column to use as labels for subset
##' @param res resolutions to cluster at
##' @param ndims number of dimensions to use for UMAP/FindNeighbors
##' @param idents clusters to subset for reclustering;
##' @param invert do you want the cluster specified by subset, or all others; default is false
##' @param nfeat number of variable features to extract
##' @return normalized and clustered seurat object
##' @author dylanmr
##' @export

quickclus <- function(seurobj, ndims = NULL, clusid = Idents(seurobj),
                      res = seq(from=0.2, to = 2, by = .2),
                      idents = NULL, invert = F,
                      nfeat = 3000) {

  # Save labels from original clustering

  seurobj[["orig_labels"]] <- Idents(seurobj)

  # Recluster

  if (is.null(idents)) {
    
    seurobj %>%
      SCTransform(variable.features.n = nfeat) %>%
      RunPCA() -> x
    
    if(is.null(ndims)) {
      ndims <- round(as.numeric(maxLikGlobalDimEst(data = x@reductions$pca[, 1:50], k = 20)))
    } else {
      ndims <- ndims
    }
  
    x %>%
      RunUMAP(dims = seq(ndims)) %>%
      FindNeighbors(dims = seq(ndims)) %>%
      FindClusters(resolution = res) -> seur
    
  } else {
    
    Idents(seurobj) <- clusid
    
    seurobj %>%
      subset(idents = idents, invert = invert) %>%
      SCTransform(variable.features.n = nfeat) %>%
      RunPCA() -> x
    
    if(is.null(ndims)) {
      ndims <- round(as.numeric(maxLikGlobalDimEst(data = x@reductions$pca[, 1:50], k = 20)))
    } else {
      ndims <- ndims
    }
    
    x %>%
      RunUMAP(dims = seq(ndims)) %>%
      FindNeighbors(dims = seq(ndims)) %>%
      FindClusters(resolution = res) -> seur
  }

  return(seur)
}
