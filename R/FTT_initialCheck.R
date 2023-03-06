#' @title Perform some initial checks
#' @description Three checks: First, the bounds of the FTT algorithm are checked. If the number of c-tuples is above the upper or below the lower bound, then the FTT algorithm does not to be used as the result is fix. Second, three properties of phylogenetic decisiveness are tested for c=4: are all triples covered? are triples with simple covered covered by different quadruples? are all tuples sufficiently covered?
#' @param data Input data created by FTT_createInput()
#' @param c Parameter indicating the size of the c-tuples. E.g. if c=4, one uploads quadruples to be tested for fixing taxon traceability. Default: 4
#' @param n Number of all taxa in the analysis
#' @param verbose Logical parameter if message should be printed, Default: F
#' @return Results of the three checks
#' @details See publication and master thesis
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{melt.data.table}},\code{\link[data.table]{copy}}
#' @rdname FTT_initialCheck
#' @export
#' @importFrom data.table melt copy
FTT_initialCheck = function(data, c=4, n, verbose = F){
  # data = InputData$data

  # Step 0: check if input is ok
  taxa = paste0("taxa",c(1:c))
  c1tuples = paste0("c1tuple",c(1:c))

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

  stopifnot(class(data$taxa1) == "integer")
  cross_tuples<-data[status=="unresolved",]
  green_tuples<-data[status=="input",]

  # Step 1: check upper bound!
  threshold1 = choose(n,c) - n + c -1

  check1 = dim(green_tuples)[1]>threshold1
  if(check1==T){
    comment1 = "Number of c-tuples above upper bound! \n     This set is always fixing-taxon traceable! You do not need to run the FTT algorithm. "
  }else{
    comment1 = "Number of c-tuples below the upper bound! \n     Please run the FTT algorithm to check for fixing-taxon traceability"
  }
  if(verbose == T){message(comment1)}

  # Step 2: check lower bound!
  threshold2 = choose(n-1,c-1)

  check2 = dim(green_tuples)[1]<threshold2
  if(check2==T){
    comment2 = "Number of c-tuples below than lower bound! \n     This set is cannot be resolved with the FTT algorithm. Please add more data!"
  }else{
    comment2 = "Number of c-tuples above the lower bound! \n     Please run the FTT algorithm to check for fixing-taxon traceability"
  }
  if(verbose == T){message(comment2)}


  # Step 3: do some more checks for c=4 (relevance for phylogenetic decisiveness)
  if((check1 == F & check2 == F & c==4)==T){
    # Step 3a: check if all c-1 tuples are covered at least once
    green_tuples_long = data.table::melt(green_tuples, id.vars = c("ctuple"), measure.vars = c1tuples)
    nr_covered_c1tuples = green_tuples_long[,.N,value]
    threshold3 = choose(n,c-1)
    check3 = dim(nr_covered_c1tuples)[1]<threshold3

    if(check3==T){
      comment3 = "There is at least on triple not covered. \n     This set is not phylogenetic decisive. \n     Please add more data!"
    }else{
      # Step 3b: check if two triples are covered by the same quadruple only
      matched = match(green_tuples_long$value,nr_covered_c1tuples$value)
      green_tuples_long[,count := nr_covered_c1tuples[matched,N]]

      green_tuples_long2 = data.table::copy(green_tuples_long)
      green_tuples_long2 = green_tuples_long2[count == 1,]
      check4 = sum(duplicated(green_tuples_long2$ctuple)) != 0

      if(check4==T){
        comment3 = "     *  all triples are covered at least once, but at least two of these triples are covered by the same quadruple. \n     This set is not phylogenetic decisive. \n     Please add more data!"

      }else{
        # Step 3c: check if tuples are covered sufficiently (see master thesis, Theorem 6)
        all_tuple_taxa<-t(combn(n,2))
        all_tuple_taxa<-as.data.frame(all_tuple_taxa)
        tuple_h<-vector(mode="numeric",length=dim(all_tuple_taxa)[1])
        for (i in 1:dim(all_tuple_taxa)[1]){
          #i=1
          tab1<-all_tuple_taxa[i,1] == green_tuples[,c(1,2,3,4)] |
            all_tuple_taxa[i,2] == green_tuples[,c(1,2,3,4)]
          dim(tab1)
          tab2<-vector(mode="numeric",length=dim(tab1)[1])
          for(j in 1:dim(tab1)[1]){
            tab2[j]<-sum(tab1[j,]==TRUE)
          }
          tuple_h[i]<-sum(tab2==2)
        }
        check5<-tuple_h<=(n-4)
        end_check5<-sum(check5==TRUE)
        if (end_check5>0){
          comment3 = "     * there is at least one tuple covered only n-3 times or less. \n     This set is not phylogenetic decisive. \n     Please add more data!"
        }else{
          comment3 = "     * all triples are covered at least once, \n     * all triples with simple coverage are covered by different c-tuples, and  \n     * all tuples are covered n-4 times. \n     Please run the FTT algorithm to check for fixing-taxon traceability."
        }

      }

    }

    comment3 = c("Checking for phylogenetic decisiveness parameters, as c=4: \n",comment3)
    if(verbose == T){message(comment3)}
  }

  if((check1 == F & check2 == F & c==4)==T){
    return(c(comment1,comment2,comment3))
  }else{
    return(c(comment1,comment2))
  }

}
