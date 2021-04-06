##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param pb_mat pseudobulk matrix 
##' @param meta metadata for adding coldata to deseq2 matrix
##' @param min_counts filter threshold: minimum number of counts detected per gene
##' @param min_samples filter threshold: minimum number of samples in which gene is detected at min_count threshold
##' @param norm_method normalize with vst or rlog; default is vst
##' @param blind When blind equals TRUE (the default), the functions will re-estimate the dispersions using only an intercept; default is false
##' 
##' @return dds object that has been filtered and normalized with vst
##' @author dylanmr
##' @export

gen_dds <- function(pb_mat, min_counts = 10, min_samples = 15, norm_method="vst", blind=F) {

  filt_glia[["celltype"]] <- Idents(filt_glia)

  filt_glia@meta.data %>%
    dplyr::distinct(cellbyid, .keep_all = T) %>% 
    dplyr::select(trt, batch, celltype, orig.ident, cellbyid) -> meta
  
  meta <- meta[match(colnames(pb_mat), meta$cellbyid),]
  rownames(meta) <- meta$cellbyid
  
  # normalize and stabilize summed counts
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pb_mat),
                                colData = meta,
                                design = ~ 0 + trt + celltype + batch)
  
  # filter genes based on defined thresholds
  keep <- Matrix::rowSums(counts(dds) >= min_counts) > min_samples
  dds <- dds[keep,]
  
  # normalize libraries
  if(norm_method == "vst") {
    norm <- vst(dds,blind = blind)
  }
  else if (norm_method == "rld") {
    norm <- rlog(dds, blind = blind)
  }
  
  return(list(dds = dds, norm = norm))
}
