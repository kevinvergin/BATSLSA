#' Phylogenetically weighted connectivity
#'
#' This analysis measures the similarity of connections for each pair of NTUs. The similarities are
#' weighted by phylogenetic distance. A profile is created for each NTU and compared, such that
#' matching, or near matching, connections score higher. It is necessary to have phylogenetic
#' distances between NTUs. This script does not provide a matrix of those values but future
#' versions of this package may include that script.
#'
#' @param inputParameter1 distmat is a 3 column matrix. Column 1 has all NTUs repeated
#' appropriately, column 2 has the 5 upstream and 5 downstream (if appropriate) NTUs paired with
#'  each corresponding NTU from column 1, and column 3 has the distance between the two
#'  NTUs from columns 1 and 2  \code{inputParameter1}
#' @param inputParameter2 boolean is a matrix of 1's and 0's from the final diagnostic scripts indicating #' which NTUs are connected (1's) \code{inputParameter2}
#'@param inputParameter3 fname is a file name for the output of the function  \code{inputParameter3}
#'@param inputParameter4 namesvec is a vector of names for the NTUs  \code{inputParameter4}
#'
#' @return output A matrix of similarity scores for the connections between each pair of NTUs weighted by their phylogenetic distance within a +/- 5 NTU window
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run script with this command:  sim_metric_csv(boolean,distmat,namesvec, fname)

# n <- length(otunames)
# n_comb <- 2*n

distadj <- function(distmat , namesvec)
{
  x <- distmat
  n <- length(namesvec)
  eta <- exp(-10)
  ymat <- matrix(nrow=n,ncol=n)
  rownames(ymat) <- namesvec
  colnames(ymat) <- namesvec

  for(i in 1:length(x[,1])) {
    ymat[as.character(x[i,1]),as.character(x[i,2])] <- x[i,3]
    ymat[as.character(x[i,1]),as.character(x[i,1])] <- 0
  }
  ymat <- exp(-5*ymat)
  ymat[is.na(ymat)] <- 0
  return(ymat)
}

 # ymatrix <- distadj()

#need a second version for combined
#this was written for comparing two datasets but is obsolete

distadj_comb <- function(xdata = x_comb, n = 2*length(otunames))
{
  eta <- exp(-10)
  ymat <- matrix(nrow=n,ncol=n)
  rownames(ymat) <- c(paste(otunames,"_000",sep=""), paste(otunames,"_200",sep=""))
  colnames(ymat) <- c(paste(otunames,"_000",sep=""), paste(otunames,"_200",sep=""))

  for(i in 1:length(x[,1])) {
    xa <- paste(as.character(x[i,1]), "_000", sep="")
    xb <- paste(as.character(x[i,2]), "_000", sep="")
    x2a <- paste(as.character(x[i,1]), "_200", sep="")
    x2b <- paste(as.character(x[i,2]), "_200", sep="")
    ymat[xa,xb] <- x[i,3]
    ymat[xb,xa] <- x[i,3]
    ymat[x2a,x2b] <- x[i,3]
    ymat[x2b,x2a] <- x[i,3]
    #ymat[xa,x2b] <- x[i,3]
    #ymat[x2b,xa] <- x[i,3]
    #ymat[x2a,xb] <- x[i,3]
    #ymat[xb,x2a] <- x[i,3]
    ymat[xa,xa] <- 0
    #ymat[xa,x2a] <- 0
    #ymat[x2a,xa] <- 0
    ymat[x2a,x2a] <- 0
  }
  ymat <- exp(-5*ymat)
  ymat[is.na(ymat)] <- 0
  return(ymat)
}

# ymatrix_comb <- distadj_comb()

#groups were not used in the final analysis so this variable can be ignored in most cases

num_sim_metric_tree <- function(boolean,group = c(),mat2 = distadj(distmat,namesvec), namesvec) {
  tmp_adjm <- boolean
  tmp_nodenames <- namesvec
  tmp_ymatrix <- mat2

  if(!is.null(group)) {
    #filter to just the group we're interested in
    tmp_adjm <- tmp_adjm[group,group]
    tmp_ymatrix <- tmp_ymatrix[group,group]
    tmp_nodenames <- tmp_nodenames[group]
  }

  #coerce nodes to numbers
  adj_nodenames <- as.character(tmp_nodenames)
  adj_nodenames <- as.numeric(adj_nodenames)

  #select those nodes which have numbered values
  tmpset <- which(!is.na(adj_nodenames))
  group_nodenames <- tmp_nodenames[tmpset]

  #smaller adjacency matrix
  tmp_adjm <- tmp_adjm[tmpset,tmpset]
  tmp_ymatrix <- tmp_ymatrix[tmpset,tmpset]

  rownames(tmp_adjm) <- group_nodenames
  colnames(tmp_adjm) <- group_nodenames

  tmp_mavem <- matrix(nrow = length(tmpset),ncol = length(tmpset))
  rownames(tmp_mavem) <- group_nodenames
  colnames(tmp_mavem) <- group_nodenames

  tmp_mavem <- tmp_adjm%*%(tmp_ymatrix)

  sim_matrix <- matrix(nrow = length(tmpset),ncol = length(tmpset))

  for(i in 1:length(tmpset)) {
    for(j in 1:length(tmpset)) {
      totpoints <- sum((tmp_mavem[i,] + tmp_mavem[j,])**2)
      diffpoints <- sum((tmp_mavem[i,] - tmp_mavem[j,])**2)

      ipoints <- sum(tmp_mavem[i,])
      jpoints <- sum(tmp_mavem[j,])
      if(ipoints == 0 | jpoints == 0) next

      sim_matrix[i,j] <- diffpoints/totpoints


    }
  }

  y <- t(1/(sim_matrix+1))

  rownames(y) <- group_nodenames
  colnames(y) <- group_nodenames

  return(y)

}

sim_metric_csv <- function(boolean,distmat,namesvec, fname) {
  tmp <- num_sim_metric_tree(boolean, namesvec)

  print(paste("Writing ", fname))
  write.csv(tmp, file=fname)
}


