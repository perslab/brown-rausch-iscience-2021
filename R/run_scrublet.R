##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author dylanmr
##' @export

run_scrublet <- function(srt, run) {

  # use reticulate to import scrublet
  scr <- import("scrublet")

  # convert srt obj sce
  sce <- Seurat::as.SingleCellExperiment(srt)

  # run scrublet with default param
  scrub <- scr$scrublet$Scrublet(t(counts(sce)))
  res <- scrub$scrub_doublets()
  names(res) <- c("dbl_score", "dbl_call")

  # check distributions of scores
  pdf(paste0(output_folder, "/scrublet/", run, "hist.pdf"))
  hist(res$dbl_score)
  dev.off()
  
  # generate dataframe with doublet results
  tibble(
    cell = colnames(sce),
    dbl_score = res$dbl_score,
    dbl_call = res$dbl_call
  ) -> dbl_results

  # return results
  return(dbl_results)

}
