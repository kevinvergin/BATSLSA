#' Checks for non-zero slope in linear models of connected NTUs
#'
#' For potential network connections between two nodal taxonomic units (NTUs),
#' verifies that the slope of a linear model is non-zero by ensuring that
#' the p value for the line is greater than 0.05
#'
#' @param inputParameter1 matrix of NTU abundances in rows and sample dates in columns \code{inputParameter1}
#' @param inputParameter2 boolean is a matrix of 1's and 0's from the LSA code indicating which NTUs are connected (1's) \code{inputParameter2}
#'@param inputParameter3 fn is a file name for the output  \code{inputParameter3}
#'
#' @return output A boolean matrix of connections that pass. Results are in the upper triangle.
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' s_matrixlm(matrix=deep_popcounts, boolean=combSPdeep,fn="lmconfirmation_200.csv")

s_matrixlm <- function(matrix, boolean,fn) {
  s_mat <- matrix(NA,nrow=nrow(boolean),ncol=ncol(boolean))
  for(i in 1:(dim(matrix)[1])) {
    for(j in 1:(dim(matrix)[1])) {
      if(i>j) next
      #only considers upper triangle
      if(i == j) next
      #ignores the diagonal

      if(boolean[i,j] < 0.5) next
      #only looks at putative correlations
      print(paste(i, "and", j))
      xval <- unlist(matrix[i,])
      yval <- unlist(matrix[j,])
      tmp <- lm(xval~yval)
      tmp2 <- anova(tmp)
      tmp3 <- tmp2$"Pr(>F)"[1]
      s_mat[i,j] <- ifelse(tmp3<=0.05,1,0)
    }
  }
  write.csv(s_mat,file=fn)
}
