#' Conservative removal of potential false positives
#'
#' Initial data screening suggested that most data did not need to be differenced seasonally first.
#' Simulations show that structure in residuals for one NTU only does not affect type I error.
#' However, structure in two NTU’s can greatly affect type I error. Check autocorrelation
#' functions for undifferenced modified data (as before). Look at row and column sums for each
#' NTU. Most NTUs have only a few acf’s with structure. Most of the acf’s with structure are
#' found in the top 5% of NTUs (sorted by the number of problem acfs). If both NTUs are in the
#' upper 5%, their type I error rate is probably elevated so they should be eliminated. Others are
#' more likely to be OK.
#' This is a two part script with commands in between to run manually. Determine an
#' appropriate cutoff (I used 5%) and then use that result in the second script. Non-automated
#' script is marked by comments below.
#'
#' @param inputParameter1 matrix of NTU abundances in rows and sample dates in columns \code{inputParameter1}
#' @param inputParameter2 boolean is a matrix of 1's and 0's from the previous linear model code indicating which NTUs are connected (1's) \code{inputParameter2}
#'@param inputParameter3 fn is a file name for the output of the first function  \code{inputParameter3}
#'@param inputParameter4 fn2 is a file name for the output of the second function  \code{inputParameter4}
#'@param inputParameter5 X is the number of connections at the cutoff  \code{inputParameter5}
#'
#' @return output A boolean matrix of connections that fail (1’s). Results are in the upper triangle. Need to subtract this matrix from the matrix of connections.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run diagnostic with this command:  s_matrixacf2(matrix=, boolean=,fn=)

deep_popcounts_trans <- read.csv("BATS_200m_timeseries_cellcounts_transposed.csv")

# use matrix=deep_popcounts_trans, boolean=lmconfirmdeep_2 in next function

s_matrixacf2 <- function(matrix, boolean,fn) {
  s_mat <- matrix(NA,nrow=nrow(boolean),ncol=ncol(boolean))
  for(i in 1:(dim(matrix)[2])) {
    for(j in 1:(dim(matrix)[2])) {
      if(i>j) next
      if(i == j) next
      if(boolean[i,j] <0.5) next
      xval <- unlist(matrix[,i])
      yval <- unlist(matrix[,j])
      tmp <- lm(xval~yval)
      tmp3 <- lm(yval~xval)
      tmp2 <- acf(tmp$resid, plot=F)
      tmp4 <- acf(tmp3$resid, plot=F)

      s_mat[i,j] <- ifelse((sum((abs(tmp2$acf[2:10]))>=0.20))>=3,1,( ifelse((sum((abs(tmp4$acf[2:10]))>=0.20))>=3,1,0)))
    }
  }
  write.csv(s_mat,file=fn)
}


#non-automated script
mat1 <- read.csv(fn)
mat1 <- mat1 [,-1] #if R adds a column of indexes
table(as.matrix(mat1))
mat1 [is.na(mat1)] <- 0 #replace NA’s with 0
rs <- rowSums(mat1)
cs <- colSums(mat1)
total1 <- rs + cs
total1 <- as.matrix(total1)
total2 <- sum(total1>0)
total3 <- total2*0.05 #can use other cutoffs, determine by graphing
sum(total1>X)  #find X by trial and error such that the sum is = total3. Fill in X in next script.
#end of non-automated script

mat2<- function(matrix, boolean,fn2,X) {
  s_mat <- matrix(NA,nrow=nrow(boolean),ncol=ncol(boolean))
  for(i in 1:(dim(boolean)[2])) {
    for(j in 1:(dim(boolean)[2])) {
      if(i>j) next
      if(i == j) next
      if(boolean[i,j] <0.5) next


      mat2[i,j] <- ifelse(matrix[i,]>=X,(ifelse(matrix[j,]>=X,1,0)),0)  #inserts 1 when two NTUs with high ACF counts are correlated. These have higher risk of being false positive.
    }
  }
  write.csv(mat2,file=fn2)
}


