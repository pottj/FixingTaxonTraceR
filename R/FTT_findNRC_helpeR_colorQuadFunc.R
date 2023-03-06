#' @title Coloring of quadruples
#' @description In this functions, the given input quadruples are colored with respect to a given coloring input. It also sums the number of colored quadruples (are there any uncolored taxa in any quadruple left?) and the sum of unique colors per quadruple (are there any quadruples with more than 3 colors?). It is a helper function for the *FTT_finNRC* algorithm.
#' @param input input data.table containing the columns taxa1, taxa2, taxa3, taxa4 and quadruple (first five columns of the data.table _data_ create by *FTT_createInput*, filtered for status == "input")
#' @param colors data.table containing all taxa with their original coding (e.g. a --> 1) and a column for the color. The four used colors are "green", "red", "blue", and "yellow". If no color is yet applied, all entries should be "nocol".
#' @return The function returns an updated input data.table. In addition the the four numeric columns, and the quadruple column, there are now the colored taxa included, plus two columns with the sums per quadruple.
#' @details needs *FTT_findNRC_helpeR_getColTaxaPerQuad* and *FTT_findNRC_helpeR_getUniColTaxaPerQuad*
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{copy}}
#' @rdname FTT_findNRC_helpeR_colorQuadFunc
#' @export
#' @importFrom data.table copy
FTT_findNRC_helpeR_colorQuadFunc = function(input,colors){
  # input = inp3
  # colors = dat1

  input2 = data.table::copy(input)
  for(i in 1:dim(colors)[1]){
    #i=1
    myTaxon = colors[i,]
    input2[taxa1 == myTaxon$NR,taxa1_col := myTaxon$color]
    input2[taxa2 == myTaxon$NR,taxa2_col := myTaxon$color]
    input2[taxa3 == myTaxon$NR,taxa3_col := myTaxon$color]
    input2[taxa4 == myTaxon$NR,taxa4_col := myTaxon$color]

  }
  input2

  # number of colors per quadruple
  uniqColPerQuad = apply(input2[,6:9], 1, FTT_findNRC_helpeR_getUniColTaxaPerQuad)

  input2[,sumUniqueColors := 0]
  input2[grepl("green",uniqColPerQuad),sumUniqueColors := sumUniqueColors + 1]
  input2[grepl("red",uniqColPerQuad),sumUniqueColors := sumUniqueColors + 1]
  input2[grepl("blue",uniqColPerQuad),sumUniqueColors := sumUniqueColors + 1]
  input2[grepl("yellow",uniqColPerQuad),sumUniqueColors := sumUniqueColors + 1]

  sumCol = apply(input2[,6:9], 1, FTT_findNRC_helpeR_getColTaxaPerQuad)

  input2[,sumColored := sumCol]
  return(input2)

}
