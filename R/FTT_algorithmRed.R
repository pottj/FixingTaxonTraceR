#' @title Alternative algorithm to test for fixing taxon traceability
#' @description Testing the input data stepwise for fixing taxa and check if the whole set can be resolved or not. The algorithm is an adaption to the *red* one, and it checks if any red c-tuple if can be solved by *c* green ones.
#' @param data Data.table as constructed by create input. All possible c-tuples given the taxon set with status information (c-tuple as input available, c-tuple not in input = unresolved). In addition, all c (c-1)-tuples possible by each c-tuple are listed.
#' @param verbose Logical parameter if message should be printed, Default: F
#' @param c Parameter indicating the size of the c-tuples. E.g. if c=4, one uploads quadruples to be tested for fixing taxon traceability. Default: 4
#' @param n Number of all taxa in the analysis
#' @return The same data.table is returned, with updated status, fixing taxon & round of resolvement.
#' @details See publication and master thesis
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{copy}}
#' @rdname FTT_algorithmRed
#' @export
#' @importFrom data.table copy
FTT_algorithmRed<-function(data, verbose = F,c=4,n){
  # data = myResults$data
  # verbose = T
  # c = 4
  # n = 6

  # Step 0: check if input is ok
  taxa = paste0("taxa",c(1:c))
  c1tuples = paste0("c1tuple",c(1:n))

  expectedNames = c(taxa,"ctuple",c1tuples,"status")
  stopifnot(names(data) %in% expectedNames)
  stopifnot(class(data$taxa1) == "integer")

  x = table(data$status)

  if(verbose == T){
    message("Using ",as.numeric(x[1])," of ",dim(data)[1]," ",c,
            "-tuples as input for algorithm (",n,
            " unique taxa). \n This leaves ",
            as.numeric(x[2])," ",c,"-tuples unsolved.")
  }

  # Step 1: Define solved and unsolved c-tuples from input data
  data2<-data.table::copy(data)
  green_ctuples<-data2[status=="input",]
  cross_ctuples<-data2[status=="unresolved",]

  # Step 2: Start repeat loop with myAlgorithm (check for fixing taxon)
  #   check each cross-c-tuples and change status if resolved
  #   repeat algorithm until all c-tuples are resolved or no changes anymore
  #   (in each round, at least one c-tuple has to be solved)
  #   (if index > dim, then there was at least one round
  #     in which no c-tuple could be solved)
  # then the set cannot be decisive

  index=0
  data2[,round := 0]
  data2[,fixingTaxa := 0]
  repeat{
    data2=FTT_algorithmRed_helpeR(data = data2,
                                  roundnumber = index,
                                  verbose = verbose,
                                  c=c,
                                  n=n)
    index=index+1
    data2

    # Loop should stop if there are no new resolved quadruples or no unresolved quadruples
    check1 = (sum(data2$round == index) == 0)
    check2 = (sum(data2$status == "unresolved") == 0)
    if (check1 | check2){
      break
    }
  }

  # Step 3: return data and result
  if(verbose == T){
    if (sum(data2$status == "unresolved") ==0){
      message("FIXING TAXON TRACEABLE")
      if(c==4)message(" It follows that the set is phylogenetically decisive")
    }else{
      x = sum(data2$status == "unresolved")
      message(paste0("NOT RESOLVABLE VIA THIS ALGORITHM! \n There are ",x," remaining cross ",c,"-tuples."))
    }
  }


  return(data2)

}
