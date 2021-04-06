##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param seur_obj
##' @param res
##' @return
##' @author dylanmr
##' @export
set_ident <- function(seur_obj, res) {

  Idents(seur_obj) <- paste0("SCT_snn_res.",res)
  return(seur_obj)

}
