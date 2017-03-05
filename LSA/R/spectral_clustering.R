#' Spectral Clustering
#'
#' This analysis measures the similarity of exact connections for each pair of NTUs. The similarities are
#' not weighted by phylogenetic distance. There are several options for spectral clustering kernels (
#' vanilladot (linear), polydot (polynomial), rbfdot (Gaussian), and others. The numbers of centers
#' should be determined using other methods, such as the pamk function in the fpc package.
#' This script takes a long time to run and could be improved by splitting tasks to separate processors.
#'
#' @param inputParameter1 matrix is the matrix of connections after the diagnostic scripts
#' are run  \code{inputParameter1}
#' @param inputParameter2 iteration is the number of bootstrap
#'  replicates to run  \code{inputParameter2}
#' @param inputParameter3 fn1 is a file name for the matrix of clustering results
#'  \code{inputParameter3}
#' @param inputParameter4 namesvec is a vector of names for the NTUs  \code{inputParameter4}
#' @param inputParameter5 centersnum is the number of clusters expected \code{inputParameter5}
#' @param inputParameter6 kerneltype is the kernel used for spectral clustering
#' \code{inputParameter6}
#' @param inputParameter7 matrix2 is the result from the bootstrapping function
#' \code{inputParameter7}
#' @param inputParameter8 fn2 is the file name for the similarity matrix for the clustering results
#' \code{inputParameter8}

#'
#' @return output 2 matrices, first the clustering results for each NTU and each iteration, and second, a #' similarity matrix based on the clustering assignments from the first part of the function
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run script with these two commands:  test1 <- bootstrap_clustering(matrix=,iteration=100,fn1="test1. #' csv",namesvec=ntunames,centersnum=7,kerneltype=rbfdot)
#' test2 <- clust_sim(matrix=, matrix2=test1,fn2="test2.csv")

no0mat <- function(matrix,namesvec){
  matrix <- as.matrix(matrix)
  rownames(matrix) <- namesvec
  colnames(matrix) <- namesvec
  startmat <- matrix[(rowSums(matrix)+colSums(matrix)>0),(rowSums(matrix)+colSums(matrix)>0)]

  return(startmat)
}



directsim <- function(startmat) {
  tmp_adjm <- startmat
  sharedconn <- tmp_adjm %*% t(tmp_adjm)
  connsbyntu <- rowSums(tmp_adjm)

  n = length(connsbyntu)
  totconn <- rep(connsbyntu,each=n) + rep(connsbyntu,n)
  dim(totconn) = c(n,n)

  #exclude
  simmat <- 2*sharedconn/totconn

  return(simmat)
}

ntulist <- function(matrix){
  ntus <- c()
  for(i in 1:nrow(matrix)){

    ifelse(sum(matrix[i,]) == 0,next, ntus <- c(ntus,rownames(matrix)[i]))
  }
  return(ntus)
}
#vector1 is from output of ntulist
samp_vector <- function(vector1){
  vector2 <- sample(vector1, round(0.9*(length(vector1))))
  return(vector2)
}

#vector1 is full ntu list, vector2 is subsample of 90% of ntus
newmatrix <- function(matrix,vector1,vector2){
  tmp_matrix <- matrix(nrow=length(vector2),ncol=length(vector2))
  rownames(tmp_matrix) <- vector2
  colnames(tmp_matrix) <- vector2
  for(i in 1:length(vector2)){
    for(j in 1:length(vector2)){
      tmp_matrix[vector2[i],vector2[j]] <- matrix[vector2[i],vector2[j]]
    }
  }
  return(tmp_matrix)
}
tmp_clust <- function(tmp_matrix,centersnum,kerneltype){
  clust <- specc(tmp_matrix,centers=centersnum,kernel="kerneltype")
  return(clust)
}
#vector1 is ntus and vector2 is sample of ntus
all_clust <- function(clust,vector1,vector2){
  clust_matrix <- matrix(nrow=1,ncol=length(vector1))
  colnames(clust_matrix) <- vector1

  for(k in 1:length(vector2)){
    clust_matrix[1,vector2[k]] <- clust@.Data[k]

  }
  return(clust_matrix)
}

#input matrix is result from bootstrap_clustering
clust_sim <- function(matrix, matrix2, fn2){
  vector1 <- ntulist(matrix)
  sim_mat <- matrix(nrow=ncol(matrix2),ncol=ncol(matrix2))
  rownames(sim_mat) <- vector1
  colnames(sim_mat) <- vector1
  for(i in 1:ncol(matrix2)){
    for(j in 1:ncol(matrix2)){
      match <- 0
      total <- 0
      for(k in 1:nrow(matrix2)){
        ifelse(is.na(matrix2[k,i]+matrix2[k,j]),next,total <- total+1)
        if((matrix2[k,i]+matrix2[k,j])/2 == (matrix2[k,i])) match <- match+1

      }
      sim_mat[i,j] <- match/total

    }
  }
  write.csv(sim_mat,file=fn2)
}

#run previous functions in proper order and iterate
bootstrap_clustering <- function(matrix,iteration,fn1,namesvec,centersnum,kerneltype){
  mat1 <- no0mat(matrix,namesvec)
  mat2 <- directsim(startmat=mat1)
  vector1 <- ntulist(mat2)
  clust_matrix_all <- matrix(nrow=iteration,ncol=length(vector1))
  colnames(clust_matrix_all) <- vector1
  for(t in 1:iteration){
    smallvector <- samp_vector(vector1)
    temp_mat <- newmatrix(matrix,vector1,vector2=smallvector)
    temp_cl <- tmp_clust(tmp_matrix=temp_mat,centersnum,kerneltype)


    for(k in 1:length(smallvector)){
      clust_matrix_all[t,smallvector[k]] <- temp_cl@.Data[k]

    }
  }
  write.csv(clust_matrix_all,file=fn1)
  return(clust_matrix_all)

}
