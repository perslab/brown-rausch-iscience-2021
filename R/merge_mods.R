##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param .x
##' @param .y
##' @param cutHeight
##' @param verbose
##' @return
##' @author dylanmr
##' @export
merge_mods <- function(dat, mods, merge_thresh) {

  merge <- mergeCloseModules(dat, dynamicColors, cutHeight = merge_thresh, verbose = 1)
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs
  moduleColors <- mergedColors
  MEs <- mergedMEs

}
