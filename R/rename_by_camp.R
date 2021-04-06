##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj
##' @return
##' @author dylanmr
##' @export
rename_by_camp <- function(seur_obj) {

  table(seur_obj$predicted.id, Idents(seur_obj)) %>% 
    as.data.frame() %>% 
    group_by(Var2) %>% 
    mutate(Percent = Freq / sum(Freq)) %>% 
    top_n(1, Percent) %>% 
    as.data.frame() %>% 
    dplyr::select(Var2, Var1) %>%
    mutate(Var2 = as.character(Var2),
           Var1 = as.character(Var1)) %>%
    deframe()

}
