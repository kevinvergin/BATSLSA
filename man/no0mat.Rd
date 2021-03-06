% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectral_clustering.R
\name{no0mat}
\alias{no0mat}
\title{Spectral Clustering}
\usage{
no0mat(matrix, namesvec)
}
\arguments{
\item{inputParameter1}{matrix is the matrix of connections after the diagnostic scripts
are run  \code{inputParameter1}}

\item{inputParameter2}{iteration is the number of bootstrap
replicates to run  \code{inputParameter2}}

\item{inputParameter3}{fn1 is a file name for the matrix of clustering results
\code{inputParameter3}}

\item{inputParameter4}{namesvec is a vector of names for the NTUs  \code{inputParameter4}}

\item{inputParameter5}{centersnum is the number of clusters expected \code{inputParameter5}}

\item{inputParameter6}{kerneltype is the kernel used for spectral clustering
\code{inputParameter6}}

\item{inputParameter7}{matrix2 is the result from the bootstrapping function
\code{inputParameter7}}

\item{inputParameter8}{fn2 is the file name for the similarity matrix for the clustering results
\code{inputParameter8}}
}
\value{
output 2 matrices, first the clustering results for each NTU and each iteration, and second, a #' similarity matrix based on the clustering assignments from the first part of the function
}
\description{
This analysis measures the similarity of exact connections for each pair of NTUs. The similarities are
not weighted by phylogenetic distance. There are several options for spectral clustering kernels (
vanilladot (linear), polydot (polynomial), rbfdot (Gaussian), and others. The numbers of centers
should be determined using other methods, such as the pamk function in the fpc package.
This script takes a long time to run and could be improved by splitting tasks to separate processors.
}
\examples{
Run script with these two commands:  test1 <- bootstrap_clustering(matrix=,iteration=100,fn1="test1. #' csv",namesvec=ntunames,centersnum=7,kerneltype=rbfdot)
test2 <- clust_sim(matrix=, matrix2=test1,fn2="test2.csv")
}
\keyword{keywords}
