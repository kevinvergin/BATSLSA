#' Spectral Cluster Sorting
#'
#' This analysis converts the results of the spectral clustering similarity matrix
#' to lists of NTUs belonging to each cluster in several steps. First, clustlist is used to convert the
#' results of the clustering analysis to a list of NTUs in each cluster (greater than 50% similarity).
#' The second step is to convert the list to a column of a matrix. A new matrix can be constructed for
#' each kernel type. The third step is to combine the cluster assignments to a new matrix using cbind (this
#' step is not included in this script). The fourth step uses clust_sim2 to make a new similarity
#' matrix based on cluster assignments. The fifth and final step is to create the final list of
#' cluster assignments (incorporating all kernel results) using clustlist2 using a greater than 75%
#' similarity cutoff.
#'
#' @param inputParameter1 matrix is the similarity matrix from the spectral clustering script
#'   \code{inputParameter1}
#' @param inputParameter2 list is the output from the clustlist script
#'   \code{inputParameter2}
#' @param inputParameter3 vector is a vector of NTU names
#'   \code{inputParameter3}
#' @param inputParameter4 matrix2 is the combined results from the clustsort script (assuming multiple kernels)
#'   \code{inputParameter4}
#' @param inputParameter5 fn3 is a file name for the clust_sim2 results
#'   \code{inputParameter5}
#' @param inputParameter6 matrix3 is the results from the clust_sim2 script
#' code{inputParameter6}
#'
#' @return output a final list of NTUs in each cluster after considering all kernels
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run scripts with these commands:  list1 <- clustlist(matrix=)
#' kernelmat <- clustsort(list=list1,vector=ntunames)
#' newmat <- clust_sim2(vector=ntunames,matrix2=cbind(kernelmat1[,2],kernelmat2[,2],etc),fn3=filename)
#' finallist <- clustlist2(matrix3=newmat)
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

clustsort <- function(list,vector){
  newvector <- c()
  for(k in 1: length(vector)){
    for(i in 1:length(list)){
      for(j in 1:length(unlist(list[i]))){
        if(vector[k]==unlist(as.numeric(list[[i]][j]))){
          newvector[k] <- i
        }

      }
    }
  }
  matrix <- cbind(vector,newvector)
  return(matrix)
}

clust_sim2 <- function(vector, matrix2, fn3){

  sim_mat <- matrix(nrow=nrow(matrix2),ncol=nrow(matrix2))
  rownames(sim_mat) <- vector
  colnames(sim_mat) <- vector
  for(i in 1:nrow(matrix2)){
    for(j in 1:nrow(matrix2)){
      match <- 0
      total <- 0
      for(k in 1:ncol(matrix2)){
        ifelse(is.na(matrix2[i,k]+matrix2[j,k]),next,total <- total+1)
        if((matrix2[i,k]+matrix2[j,k])/2 == (matrix2[i,k])) match <- match+1

      }
      sim_mat[i,j] <- match/total

    }
  }
  write.csv(sim_mat,file=fn3)
}



clustlist2 <- function(matrix3,i=1,n=1,clusts=list()){

  ntus <- c()
  while(i <= ncol(matrix3)){
    g <- 0
    # need to create the first/new list, length starts at 0
    if(length(unlist(clusts[n]))==0){
      for(j in 1:ncol(matrix3)){
        if(i>=j) next
        if(length(unlist(clusts[n]))==0){
          if(matrix3[i,j] > 0.75){
            clusts[n] <- colnames(matrix3)[i]
            clusts[[n]] <- c(clusts[[n]],colnames(matrix3)[j])
            ntus <- c(ntus,colnames(matrix3)[i])
            ntus <- c(ntus,colnames(matrix3)[j])
            if(length(unlist(clusts[n])) >0) next
          }
        }
        # if the first/new list is started, can add to it
        if(length(unlist(clusts[n]))>=2){

          if(matrix3[i,j] > 0.75){
            clusts[[n]] <- c(clusts[[n]],colnames(matrix3)[j])
            ntus <- c(ntus,colnames(matrix3)[j])
          }
        }
      }
    }
    # if the first list is still empty after the first row i, go to next i. g is a dummy variable to determine loops to run. If first list is created already, need to check next i against previous entries and skip if already listed.

    if(length(unlist(clusts[n]))==0){
      if(length(ntus)>0){
        for(m in 1:length(ntus)){
          if(i < ncol(matrix3)){


            if(as.numeric(colnames(matrix3))[i+1] == sort(as.numeric(ntus))[m]){
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
        if(i < ncol(matrix3)){
          if(as.numeric(colnames(matrix3))[i+1] == sort(as.numeric(ntus))[m]){
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
        if(i < ncol(matrix3)){
          i <- i+1
          n <- n+1
          g <- 0
        }
      }




    }

  }
  return(clusts)
}


