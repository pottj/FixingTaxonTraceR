#' @title Helper function for the *red* algorithm for fixing taxon traceability.
#' @description Split data to input and unresolved c-tuples and try for each unresolved one if a fixing taxon can be found.
#' @param data Data.table as constructed by create input + count and fixing taxon. All possible c-tuples given the taxon set with status information (c-tuple as input available, c-tuple not in input = unresolved, c-tuples already solved = resolved). In addition, all c (c-1)-tuples possible by each c-tuple are listed. For resolved c-tuples, used fixing taxon and round of resolvement is listed.
#' @param roundnumber Index for round
#' @param verbose Logical parameter if message should be printed, Default: F
#' @param c Parameter indicating the size of the c-tuples. E.g. if c=4, one uploads quadruples to be tested for fixing taxon traceability. Default: 4
#' @param n Number of all taxa in the analysis
#' @return The same data.table is returned, with updated status, fixing taxon & round
#' @details See publication and master thesis
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[foreach]{foreach}}
#'  \code{\link[data.table]{data.table-package}},\code{\link[data.table]{copy}},\code{\link[data.table]{rbindlist}}
#' @rdname FTT_algorithmRed_helpeR
#' @export
#' @importFrom foreach foreach
#' @importFrom data.table data.table copy rbindlist
FTT_algorithmRed_helpeR<-function(data,roundnumber,verbose = F,c=4,n){
  # data = data.table::copy(data2)
  # roundnumber = index
  # verbose = T
  # c = 4
  # n = 6

  # Step 0: check if input is ok
  taxa = paste0("taxa",c(1:c))
  c1tuples = paste0("c1tuple",c(1:n))

  expectedNames = c(taxa,"ctuple",c1tuples,"status","round","fixingTaxa")
  stopifnot(names(data) %in% expectedNames)
  stopifnot(class(data$taxa1) == "integer")
  round_old = max(data$round)

  # Step 1: split quadruples in resolved and unresolved ones
  data_solved = data[status != "unresolved"]
  data_unresolved = data[status == "unresolved"]

  # Step 2: check every unresolved c-tuple for fixing taxon
  myTab = foreach::foreach(i = 1:dim(data_unresolved)[1]) %do%{
    # i=1

    # Step 2.1: get all (c-1)-tuple for unresolved c-tuple
    filt = grepl("c1tuple",names(data_unresolved))
    triples_CQ<-data_unresolved[i,filt,with=F]
    triples_CQ = unique(unlist(triples_CQ))

    # Step 2.2: check overlap with resolved (c-1)-tuple
    dumTab = foreach::foreach(j = 1:c)%do%{
      # j=1
      vgl = is.element(data_solved[,get(paste0("c1tuple",j))],triples_CQ)
      vgl
    }

    filt = dumTab[[1]]
    for(j in 2:c){
      #j=2
      vgl = dumTab[[j]]
      filt<-filt | vgl
    }

    # Step 2.3: get resolved c-tuples that contain at least one of the four tested (c-1)-tuple
    data_pos<-data_solved[filt,]

    if(dim(data_pos)[1]==0){
      # Step 2.5: return best taxa, taxa count & unresolved c-tuple
      res = data.table::data.table(unres_tuples = data_unresolved[i,ctuple],
                                   best_fixTaxa = 0,
                                   count = 0)
    }else{
      # Step 2.4: search for possible fixing taxa
      # using only taxa in data_pos, which are not in the unresolved c-tuple
      # Count the most common taxa only!
      taxa_pos2 = data.table::copy(data_pos)
      taxa_pos2 = taxa_pos2[,1:4,with=F]
      taxa_pos3 = unlist(taxa_pos2)
      taxa_pos3 = as.numeric(taxa_pos3)

      data_unresolved2 = data.table::copy(data_unresolved)
      data_unresolved2 = data_unresolved2[i,1:4,with=F]
      data_unresolved3 = unlist(data_unresolved2)
      data_unresolved3 = as.numeric(data_unresolved3)

      taxa_pos2 = taxa_pos3[!is.element(taxa_pos3,data_unresolved3)]
      tab = data.table::data.table(taxa_ID = taxa_pos2)
      tab2 = tab[,.N, by=taxa_ID]
      tab3 = tab2[N == max(N),]
      if(dim(tab3)[1]>1){tab3 = tab3[1,]}

      # Step 2.5: return best taxa, taxa count & unresolved c-tuple
      res = data.table::data.table(unres_tuples = data_unresolved[i,ctuple],
                                   best_fixTaxa = tab3$taxa_ID,
                                   count = tab3$N)

    }

    res
  }
  myTab = data.table::rbindlist(myTab)
  myTab

  # Step 3: if count of taxa >=c it can be considered a fixing taxon, and the c-tuple is resolved
  stopifnot(myTab$unres_tuples == data_unresolved$ctuple)
  filt = myTab$count>=c
  data_unresolved[filt,status := "resolved"]
  data_unresolved[filt,fixingTaxa := myTab[filt,best_fixTaxa]]
  data_unresolved[filt,round := roundnumber + 1]
  data_unresolved

  # Step 4: return data
  data3 = rbind(data_solved,data_unresolved)
  table(data3$status)
  tab4<-data3[,.N,by=status]
  n_unresolved_new = tab4[status == "unresolved",N]
  if (length(n_unresolved_new)==0){n_unresolved_new = 0}
  n_unresolved_old = dim(data_unresolved)[1]
  n_diff = n_unresolved_old - n_unresolved_new
  round_new = max(data3$round)
  if(round_new == round_old){round_new = round_new + 1}

  if(verbose == T){message("In round #",round_new,", ",n_diff," ",c,"-tuples could be resolved ...")}
  return(data3)
}
