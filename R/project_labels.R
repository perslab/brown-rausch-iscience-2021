##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param ref
##' @param query
##' @return
##' @author dylanmr
##' @export

project_labels <- function(ref, query) {

    anch <- FindTransferAnchors(reference = ref, query = query, normalization.method = "SCT", dims = 1:30)
    predictions <- TransferData(anchorset = anch, refdata = Idents(ref))
    query <- AddMetaData(query, metadata = predictions)
    
    return(query)

}
