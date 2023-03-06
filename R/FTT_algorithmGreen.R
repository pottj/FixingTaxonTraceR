#' @title Algorithm to test for fixing taxon traceability
#' @description Testing the input data stepwise for fixing taxa and check if the whole set can be resolved or not. The algorithm follows strictly the Mathematica template from Mareike Fischer by checking each green c-tuple if it can be used to solve a red one. There are two counters: first a simple counter adding up all solved c-tuples (green + newGreens). If this reaches the maximal number of c-tuples, the while loops stops. Second, I use a count to track the number of while and for loops used.
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
#'  \code{\link[foreach]{foreach}}
#' @rdname FTT_algorithmGreen
#' @export
#' @importFrom data.table copy
#' @importFrom foreach foreach
FTT_algorithmGreen<-function(data, verbose = F,c=4,n){
  # data = example_3.4$data
  # verbose = T

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

  # Step 1: define parameters
  # X = unique(unlist(data2[,c("taxa1","taxa2","taxa3","taxa4")]))
  # n = length(X)
  data2<-data.table::copy(data)
  data2[,round := 0]
  data2[,counter := 0]
  data2[,fixingTaxa := 0]
  data2[status != "input" ,round := NA]
  data2[status != "input" ,counter := NA]
  data2[status != "input" ,fixingTaxa := NA]

  X = c(1:n)
  maxCounter = choose(n,c)
  newGreens = data2[status =="input",]
  allGreens = data2[status =="input",]
  counter = dim(newGreens)[1]
  count = 0
  data_unresolved = data2[status !="input",]

  # Step 2: While loop
  while(counter < maxCounter & dim(newGreens)[1]!=0){
    quad = newGreens[1,]
    filt = grepl("taxa",names(quad))
    inputTaxa  = quad[,filt,with=F]
    inputTaxa = unlist(inputTaxa)
    X2 = X[!is.element(X,inputTaxa)]

    myTab = foreach::foreach(i = 1:length(X2))%do%{
      # i=1
      x = X2[i]
      test_quad1 = c(x,inputTaxa[-1])
      test_quad1 = test_quad1[order(test_quad1)]
      test_quad1 = paste(test_quad1,collapse="_")
      test_quad2 = c(x,inputTaxa[-2])
      test_quad2 = test_quad2[order(test_quad2)]
      test_quad2 = paste(test_quad2,collapse="_")
      test_quad3 = c(x,inputTaxa[-3])
      test_quad3 = test_quad3[order(test_quad3)]
      test_quad3 = paste(test_quad3,collapse="_")
      test_quad4 = c(x,inputTaxa[-4])
      test_quad4 = test_quad4[order(test_quad4)]
      test_quad4 = paste(test_quad4,collapse = "_")

      data_pos = allGreens[is.element(ctuple,c(test_quad1,test_quad2,test_quad3,test_quad4))]

      if(dim(data_pos)[1]==3){
        filt = is.element(data_unresolved$ctuple,c(test_quad1,test_quad2,test_quad3,test_quad4))
        if(data_unresolved[filt,status] == "unresolved"){
          data_unresolved[filt,status := "resolved"]
          data_unresolved[filt,round := count + 1]
          data_unresolved[filt,counter := counter + 1]
          filt3 = grepl("taxa",names(data_unresolved))
          resolvedTaxa = data_unresolved[filt,filt3,with=F]
          resolvedTaxa = unlist(resolvedTaxa)
          FT = inputTaxa[!is.element(inputTaxa,resolvedTaxa)]
          data_unresolved[filt,fixingTaxa := FT]
          newGreens = rbind(newGreens,data_unresolved[filt,])
          allGreens = rbind(allGreens,data_unresolved[filt,])
          data_unresolved = data_unresolved[!filt,]
          counter = counter + 1
        }
      }
      count = count + 1
    }

    # get new green and new count
    newGreens = newGreens[!is.element(ctuple,quad),]


  }
  count
  counter


  # Step 3: Check result
  if(verbose == T){
    message("It took ",count," steps to come to a conclusion ...")
  }

  if(verbose == T){
    if (counter == maxCounter){
      message("FIXING TAXON TRACEABLE")
      if(c==4)message(" It follows that the set is phylogenetically decisive")
    }else{
      x = sum(data_unresolved$status == "unresolved")
      message(paste0("NOT RESOLVABLE VIA THIS ALGORITHM! \n There are ",x," remaining cross ",c,"-tuples."))
    }
  }

  finData = rbind(allGreens,data_unresolved)
  finData[maxCounter,round := count]
  return(finData)

}
