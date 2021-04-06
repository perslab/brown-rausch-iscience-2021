##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param srt
##' @return
##' @author dylanmr
##' @export
gen_model_obj <- function(srt, group) {

  datlist <- SplitObject(srt,split.by = group)
  
  as_tibble(names(datlist)) %>%
    mutate(data = map(datlist, ~ norm_by_ct(.x)),
           meta = map(datlist, ~ .x[[]]),
           meta = map(meta, ~ dplyr::mutate(.x, trt = dplyr::if_else(is.na(trt), true = "R50", false = trt)))) -> test
  
  return(test)

}
