##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param path
##' @return
##' @author dylanmr
##' @export

gen_campbell_ref <- function(path = NULL) {

  path = "/projects/dylan/rausch_kapel_2020/data/campbell_all_clust_full.RDS.gz"
  
  x <- readRDS(path)
  
  filt_idents <- as.character(unique(Idents(x))[grepl("^n.*", unique(Idents(x)))])
  
  x %>%
    subset(idents = filt_idents, invert=T) %>% 
    SCTransform() %>%
    RunPCA() -> glia
  
  x %>%
    subset(idents = filt_idents) %>% 
    SCTransform() %>%
    RunPCA() -> neur
  
  return(list(glia = glia, neur=neur))
}
