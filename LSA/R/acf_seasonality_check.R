#' Checks non-normal distribution of residuals from linear model by examining the autocorrelation
#' function
#'
#' Need a dataset with no gaps in the time series. The autocorrelation function is sensitive to the
#' order of the variables in the linear model so generates data for both cases. Difference
#' the data (lag=12 for monthly data) to pre-whiten the data.
#'
#' @param inputParameter1 matrix of NTU abundances in rows and sample dates in columns \code{inputParameter1}
#' @param inputParameter2 boolean is a matrix of 1's and 0's from the previous linear model code indicating which NTUs are connected (1's) \code{inputParameter2}
#'@param inputParameter3 fn is a file name for the output  \code{inputParameter3}
#'
#' @return output A boolean matrix of connections that fail (1’s). Results are in the upper triangle. Need to subtract this matrix from the matrix of connections.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' Run diagnostic with this command:  s_matrixacf(matrix=, boolean=,fn=)

s_matrixacf <- function(matrix, boolean,fn) {
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
