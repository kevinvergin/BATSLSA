#' Spectral Cluster Sorting
#'
#' This analysis converts the results of the spectral clustering similarity matrix
#' to lists of NTUs belonging to each cluster.
#'
#' @param inputParameter1 matrix is the similarity matrix from the spectral clustering script
#'   \code{inputParameter1}


#'
#' @return output a list of NTUs in each cluster (arbitrarily numbered)
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run script with this command:  list1 <- clustlist(matrix=)
#'

clustlist <- function(matrix,i=1,n=1,clusts=list()){

  ntus <- c()
  while(i <= ncol(matrix)){
    g <- 0
    # need to create the first/new list, length starts at 0
    if(length(unlist(clusts[n]))==0){
      for(j in 1:ncol(matrix)){
        if(i>=j) next
        if(length(unlist(clusts[n]))==0){
          if(matrix[i,j] > 0.5){
            clusts[n] <- colnames(matrix)[i]
            clusts[[n]] <- c(clusts[[n]],colnames(matrix)[j])
            ntus <- c(ntus,colnames(matrix)[i])
            ntus <- c(ntus,colnames(matrix)[j])
            if(length(unlist(clusts[n])) >0) next
          }
        }
        # if the first/new list is started, can add to it
        if(length(unlist(clusts[n]))>=2){

          if(matrix[i,j] > 0.5){
            clusts[[n]] <- c(clusts[[n]],colnames(matrix)[j])
            ntus <- c(ntus,colnames(matrix)[j])
          }
        }
      }
    }
    # if the first list is still empty after the first row i, go to next i. g is a dummy variable to determine loops to run. If first list is created already, need to check next i against previous entries and skip if already listed.

    if(length(unlist(clusts[n]))==0){
      if(length(ntus)>0){
        for(m in 1:length(ntus)){
          if(i < ncol(matrix)){


            if(as.numeric(colnames(matrix))[i+1] == sort(as.numeric(ntus))[m]){
              i <- i+1
              g <- 3
            }



          }
        }
      } else {
        i <- i+1
        g <- 3
      }
      i <- i+1
    }
    #check next i against previous entries on all lists and skip if already listed
    if(length(unlist(clusts[n]))>0){
      for(m in 1:length(ntus)){
        if(i < ncol(matrix)){
          if(as.numeric(colnames(matrix))[i+1] == sort(as.numeric(ntus))[m]){
            #if(i < ncol(matrix)){
            i <- i+1
            g <- 1
          }
        }


      }

      if(g==1){
        i <- i+1
        n <- n+1
        g <- 2
      }
      if(g == 0){
        if(i < ncol(matrix)){
          i <- i+1
          n <- n+1
          g <- 0
        }
      }




    }

  }
  return(clusts)
}
