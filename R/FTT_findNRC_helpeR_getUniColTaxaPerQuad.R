#' @title Get number of unique colors per quadruple
#' @description Helper function for *FTT_findNRC_helpeR_colorQuadFunc* to sum up all unique colors of a quadruple
#' @param x vector with four character elements, either "green", "red", "blue", "yellow" or "nocol".
#' @return character containing all colors, seperated by "|"
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname FTT_findNRC_helpeR_getUniColTaxaPerQuad
#' @export
FTT_findNRC_helpeR_getUniColTaxaPerQuad = function(x){
  y = unique(x)
  y = paste(y,collapse = "|")
  return(y)
}
