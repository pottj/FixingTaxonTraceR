#' @title Get number of colored taxa per quadruple
#' @description Helper function for *FTT_findNRC_helpeR_colorQuadFunc* to sum up all entries per row that are colored
#' @param x vector with four character elements, either "green", "red", "blue", "yellow" or "nocol".
#' @return Sum of colored taxa for the input quadruple
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FTT_findNRC_helpeR_getColTaxaPerQuad
#' @export
FTT_findNRC_helpeR_getColTaxaPerQuad = function(x){
  y = sum(x != "nocol" )
  return(y)
}
