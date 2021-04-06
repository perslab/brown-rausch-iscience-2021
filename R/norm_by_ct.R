##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param srt seurat object
##' @param clusid how you want clusters split

##' @return
##' @author dylanmr
##' @export

norm_by_ct <- function(srt, clusid) {

  SCTransform(srt, verbose=F) %>%
    GetAssayData(., slot = "scale.data") %>%
    t() %>%
    as_tibble(rownames = NA) %>%
    janitor::clean_names()

}
