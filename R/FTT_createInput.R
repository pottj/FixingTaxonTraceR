#' @title Create Input for FixingTaxonTraceR Package
#' @description The function takes a text file with taxon lists (c-tuples or more) and transforms them to target format (only c-tuples, numbered, ordered)
#' @param fn Path to text file containing the taxon list. Taxa within can be numbers or characters. If the taxon lists have different length, missing entries will be filled with NA (if input is numeric) or "" (if input are characters).
#' @param sepSym Character or symbol used to separate the taxa, e.g. "_" or ","
#' @param c Parameter indicating the size of the c-tuples. E.g. if c=4, one uploads quadruples to be tested for fixing taxon traceability. Default: 4
#' @return The output is a list containing the 5 data sets:
#'
#' - input_raw: data as given in the text file.
#' - input_ctuples: interim result, data transformed into c-tuples only
#' - input_ordered: data transformed to c-tuples, and taxa ID forced to numeric and ordered hierarchically
#' - data: all possible c-tuples given the taxon set with status information (c-tuple as input available, c-tuple not in input = unresolved). In addition, all possible (c-1)-tuples by each c-tuple are listed
#' - taxa: data table used for transformation, taxaID denotes the original input taxaID (as in input_raw & input_quadruples), NR is the ordered number of this taxon (as in input_ordered & data)
#'
#' @details
#'
#' 1) Get number of taxon lists (trees), unique taxa count, and transformation matrix
#' 2) Split larger list into quadruples
#' 3) Force numeric & ordered taxa (usind transformation matrix)
#' 4) Get all possible quadruples given the number of unique input taxa
#' 5) Get status information & return list
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{setDTthreads}},\code{\link[data.table]{fread}},\code{\link[data.table]{data.table-package}},\code{\link[data.table]{setorder}},\code{\link[data.table]{as.data.table}},\code{\link[data.table]{rbindlist}},\code{\link[data.table]{setattr}},\code{\link[data.table]{copy}}
#'  \code{\link[foreach]{foreach}}
#' @rdname FTT_createInput
#' @export
#' @importFrom data.table setDTthreads fread data.table setorder as.data.table rbindlist setnames copy
#' @importFrom foreach foreach
FTT_createInput<-function(fn, sepSym, c=4){
  # fn = "inst/extdata/example_2_Decisive.txt"
  # sepSym = ","
  # c=4

  # Step 0: load data
  data.table::setDTthreads(1)
  input<-data.table::fread(fn,sep=sepSym, fill=T,header=F)

  # Step 1: count columns & rows
  nCol = dim(input)[2]
  nRow = dim(input)[1]
  UniqueTaxa = unique(unlist(input))
  UniqueTaxa = UniqueTaxa[!is.na(UniqueTaxa)]
  UniqueTaxa = UniqueTaxa[UniqueTaxa != ""]
  myTrafoMatrix = data.table::data.table(taxaID = UniqueTaxa)
  names(myTrafoMatrix) = "taxaID"
  data.table::setorder(myTrafoMatrix,taxaID)
  myTrafoMatrix[,NR:=1:dim(myTrafoMatrix)[1]]

  stopifnot(nCol>=c)
  message("Input contains ",nRow," sets with ",dim(myTrafoMatrix)[1],
          " different taxa. \nThe largest set has ",nCol," taxa.")

  # Step 2: get c-tupel format
  # if there are more than c taxa, I can assume that all
  # possible c-tuples are known.
  # All those c-tuples are hence added to the input data

  myData = foreach::foreach(j = 1:nRow)%do%{
    # j=1
    myRow = input[j,]

    UniqueTaxa2= unique(unlist(myRow))
    UniqueTaxa2 = UniqueTaxa2[!is.na(UniqueTaxa2)]
    UniqueTaxa2 = UniqueTaxa2[UniqueTaxa2 != ""]
    myTrafoMatrix2 = data.table::data.table(taxaID = UniqueTaxa2)
    names(myTrafoMatrix2) = "taxaID"
    data.table::setorder(myTrafoMatrix2,taxaID)
    myTrafoMatrix2[,NR:=1:dim(myTrafoMatrix2)[1]]

    all_ctuples<-t(combn(dim(myTrafoMatrix2)[1],c))
    all_ctuples<-data.table::as.data.table(all_ctuples)
    names(all_ctuples) = paste0("taxa",c(1:c))

    all_ctuples2 = all_ctuples

    for(k in dim(myTrafoMatrix2)[1]:1){
      # k=4
      filt = all_ctuples2 == myTrafoMatrix2[k,NR]
      all_ctuples2[filt] = myTrafoMatrix2[k,taxaID]
    }
    ctuple = c()
    for(l in 1:dim(all_ctuples2)[1]){
      #l=1
      myRow2 = all_ctuples2[l]
      ctuple_dummy = paste(myRow2, collapse ="_")
      ctuple = c(ctuple,ctuple_dummy)
    }
    all_ctuples2[,ctuple:=ctuple]
    all_ctuples2
  }
  myData = data.table::rbindlist(myData)
  myData = myData[!duplicated(ctuple),]

  # Step 3: force numeric & ordered taxa
  myData2 = data.table::copy(myData)
  for(k in dim(myTrafoMatrix)[1]:1){
    # k=4
    filt = myData2 == myTrafoMatrix[k,taxaID]
    myData2[filt] = myTrafoMatrix[k,NR]
  }
  myData2

  # Step 4: get all possible quadruples given the number of input taxa
  n = dim(myTrafoMatrix)[1]
  all_ctuples<-t(combn(n,c))
  all_ctuples<-data.table::as.data.table(all_ctuples)
  names(all_ctuples) = paste0("taxa",c(1:c))

  # Step 5: get overview table with all taxa, triple and quadruples
  myTestTab = foreach::foreach(j = 1:dim(all_ctuples)[1])%do%{
    # j=1
    myRow = all_ctuples[j,]

    # get ctuple
    ctuple = paste0(myRow, collapse = "_")
    myRow[,ctuple := ctuple]

    # get (c-1)tuple
    myTestTab2 = foreach::foreach(k = c:1)%do%{
      # k=1
      myRow2 = data.table::copy(myRow)
      myRow2 = myRow2[,1:c,with=F]
      myRow2 = myRow2[,-k,with=F]
      ctuple_1 = paste0(myRow2, collapse = "_")
      myRow[,triple := ctuple_1]
      data.table::setnames(myRow,"triple",paste0("c1tuple",k))
      myRow
    }
    myRow

  }
  all_ctuples = data.table::rbindlist(myTestTab)

  data_all<-data.table::copy(all_ctuples)

  filt = is.element(data_all$ctuple,myData2$ctuple)
  table(filt)
  data_all[filt,status:="input"]
  data_all[!filt,status:="unresolved"]

  # return input as list
  myData[,ctuple := NULL]
  myResults<-list(input_raw = input,
                  input_ctuples = myData,
                  input_ordered = myData2,
                  data = data_all,
                  taxa = myTrafoMatrix)
  return(myResults)
}
