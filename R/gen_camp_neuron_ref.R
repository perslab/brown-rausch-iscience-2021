##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param path
##' @return
##' @author dylanmr
##' @export

gen_camp_neuron_ref <- function(path = "/projects/dylan/rausch_kapel_2020/data/campbell_all_clust_full.RDS.gz") {

  x <- readRDS(path)
  
  filt_idents <- as.character(unique(Idents(x))[grepl("^n.*", unique(Idents(x)))])
  
  x %>%
    subset(idents = filt_idents) %>% 
    SCTransform() %>%
    RunPCA() -> neur

  neur[["celltype"]] <- Idents(neur)
  neur[["celltype"]] <- sapply(strsplit(as.character(neur[["celltype"]][,1]), "[.]"), "[",2)
  Idents(neur) <- "celltype"
    
  return(neur)
}