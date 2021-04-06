##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seurobj seurat object
##' @param features list of features to build scores from
##' @param names matched vector of names for each set of features
##' @param group what group do you want to classify (should be stored in metadata column of seurat object)
##' @param k number of expected gaussian distributions
##' @return
##' @author dylanmr
##' @export

library(mixtools)

classify_by_geneset <- function(seurobj, features, 
                                names, group, k=2) {
  
  # set groups to classify
  Idents(seurobj) <- group
  
  # score each geneset
  seurobj <- AddModuleScore(seurobj, features = features, name = names)
  
  # update name to access scores 
  acc_name <- paste0(names,1)
  
  # calculate average score across all groups
  tibble(
    score = unlist(seurobj[[acc_name]]), 
    cluster = Idents(seurobj)
    ) %>%
    group_by(cluster) %>% 
    summarise(mean=mean(score)) -> score

  # run gaussian micture model to identify k distributions
  mixmdl <- normalmixEM(score$mean, k = k)
  
  # classify each group according to posterior values
  data.frame(id = levels(score$cluster), mixmdl$posterior) %>%
    mutate(assignment = if_else(comp.1 > comp.2, true = "decr", false= "incr")) -> results
  
  # match grouping id to outcome and return
  results %>%
    left_join(seurobj[[]], by= c("id" = group)) %>%
    distinct(id, .keep_all = T) %>%
    dplyr::select(id, assignment)
  
}
