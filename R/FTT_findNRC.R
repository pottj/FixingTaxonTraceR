#' @title Find no-rainbow-coloring for input quadruples
#' @description This function searches for a no-rainbow 4-coloring of the input quadruples. If there is one, then it represents one 4-way partition that is not covered by the input quadruples - hence the qaudruple set is not phylogenetic decisive.
#'
#' **IMPORTANT NOTE**: This function is only meaningful for c=4.
#' @param data The input data is the list created by *FTT_createInput*, both the data table *taxa* and *data* is used
#' @param verbose Logical parameter if all messages should be printed, Default: F
#' @param verbose2 Logical parameter if final message should be printed, Default: T
#' @return A list containing three sets:
#' - the coloring per taxa
#' - the colored quadruples
#' - the result: _fail_ indicates decisiveness (no coloring was a no-rainbow 4-coloring), while _NRC_ indicates the first no-rainbow 4-coloring and the reason to call the set not-decisive
#' @details This algorithm was adapted from Algorithm 2 (No-Rainbow 4-Coloring of H) of [Ghazaleh Parvini, Katherine Braught, and David Fernandez-Baca](https://ieeexplore.ieee.org/document/9616390) (DOI: 10.1109/TCBB.2021.3128381)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[data.table]{copy}}
#' @rdname FTT_findNRC
#' @export
#' @importFrom data.table copy
FTT_findNRC<-function(data, verbose = F, verbose2 = T){
  # data = input
  dat1 = data$taxa
  inp1 = data$data[status=="input",1:5]

  n = dim(dat1)[1]
  myI = floor(n/4)
  dat1[,color := "nocol"]

  inp2 = FTT_findNRC_helpeR_colorQuadFunc(input = inp1,colors = dat1)

  counter = 0

  for(i in 1:myI){
    #i=1

    # step1: uncul
    x = c(1:n)
    dat1[,color := "nocol"]

    a = t(combn(n,i))
    a1 = dim(a)[1]

    for(k in 1:a1){
      #k=1
      a2 = a[k,]
      x0 = x[a2]
      x1 = x[!is.element(x,x0)]

      dat2 = data.table::copy(dat1)
      dat2[NR %in% x0, color := "green"]

      myJ = floor((n-i)/3)

      for(j in 1:myJ){
        #j=1
        b = t(combn(length(x1),j))
        b1 = dim(b)[1]

        for(l in 1:b1){
          #l=1
          b2 = b[l,]
          x2 = x1[b2]
          x3 = x1[!is.element(x1,x2)]

          dat3 = data.table::copy(dat2)
          dat3[NR %in% x2, color := "red"]

          # coloring quadruples
          if(verbose==T) message("Testing i=",i,", k=",k,", j=",j,", and l=",l," ...")
          inp3 = FTT_findNRC_helpeR_colorQuadFunc(input = inp2,colors = dat3)

          if(max(inp3$sumColored)<2){
            x4 = x3[1]
            x5 = x3[!is.element(x3,x4)]
            dat3[NR==x4, color := "blue"]
            dat3[NR %in% x5, color := "yellow"]
            inp3 = FTT_findNRC_helpeR_colorQuadFunc(input = inp2,colors = dat3)
            myRes = list(coloring = dat3,
                         quadruples_colored = inp3,
                         result = "NRC")
            #print(myRes$result)

          }else{
            inp4 = data.table::copy(inp3)
            inp4 = inp4[sumColored==2 & sumUniqueColors==2,]
            inp4 = inp4[1,]
            myFilt = grep("nocol",inp4[1,])
            myFilt = myFilt - 5

            x4 = inp4[1,myFilt,with=F]
            dat3[NR %in% x4, color := "blue"]
            inp3 = FTT_findNRC_helpeR_colorQuadFunc(input = inp3,colors = dat3)

            filt = inp3$sumUniqueColors==3 & inp3$sumColored==3
            filt2 = sum(filt)
            while(filt2>0){
              inp4 = data.table::copy(inp3)
              inp4 = inp4[sumColored==3 & sumUniqueColors==3,]
              inp4 = inp4[1,]
              myFilt = grep("nocol",inp4[1,])
              myFilt = myFilt - 5

              x4 = inp4[1,myFilt,with=F]
              dat3[NR %in% x4, color := "blue"]
              inp3 = FTT_findNRC_helpeR_colorQuadFunc(input = inp3,colors = dat3)

              filt = inp3$sumUniqueColors==3 & inp3$sumColored==3
              filt2 = sum(filt)
            }

            if(sum(inp3$sumColored==4)==dim(inp3)[1]){
              myRes = list(coloring = dat3,
                           quadruples_colored = inp3,
                           result = "fail")
              #print(myRes$result)
            }else{
              dat3[color == "nocol", color := "yellow"]
              inp3 = FTT_findNRC_helpeR_colorQuadFunc(input = inp3,colors = dat3)
              myRes = list(coloring = dat3,
                           quadruples_colored = inp3,
                           result = "NRC")
              #print(myRes$result)
            }

          }
          if(i==1 & k==1 & j==1 & l==1) counter = counter +1
          if(verbose==T) message("        this coloring resulted in: ",myRes$result,"\n        counter at ",counter)
          finalFilt = myRes$result == "fail" & counter<2
          #print(finalFilt)
          #if(finalFilt == TRUE) rm(myRes)
          if(finalFilt == FALSE) break

        }
        if(finalFilt == FALSE) break
      }
      if(finalFilt == FALSE) break
    }
    if(finalFilt == FALSE) break
  }

  # Step 3: return data and result
  if(verbose2 == T){
    if (myRes$result =="fail"){
      message("THERE IS ONLY RAINBOW COLORING POSSIBLE OR 3-WAY COLORING")
      message(" It follows that the set is phylogenetically decisive")
    }else{
      message(paste0("THERE IS AT LEAST ONE NO RAINBOW COLORING IS POSSIBLE"))
      message(" It follows that the set is not phylogenetically decisive")
    }
  }
  return(myRes)

}
