extractmatrices <- function(LCO, delayrange=c(1), fdr=TRUE) {
  corn_d = LCO$corn_d
  corns_d = LCO$corns_d

  if(fdr) {
    #adjust p-value according to Benjamini-Hochberg to control FDR
    dim(corns_d) = c()
    corns_d = p.adjust(corns_d,method="BH",n=sum(!is.na(corns_d)))
    dim(corns_d) = dim(LCO$corns_d)
  }
  otu1 = dim(corns_d)[[1]]
  otu2 = dim(corns_d)[[2]]
  delay = dim(corns_d)[[3]]
  finalmatrix = matrix(0,nrow=otu1,ncol=otu2)
  delaymatrix = matrix(NA,nrow=otu1,ncol=otu2)
  cormatrix = matrix(0,nrow=otu1,ncol=otu2)
  for(i in 1:otu1) {
    for(j in 1:otu2) {
      smallestp = which.min(corns_d[i,j,])
      effectivp = which.min(corns_d[i,j,delayrange])
      if(length(effectivp)==0) next
      if(corns_d[i,j,effectivp]<=0.01) {
        delaymatrix[i,j] = smallestp-1
        cormatrix[i,j] = corn_d[i,j,effectivp]
      }
    }
  }
  signmatrix = sign(cormatrix)
  finalmatrix = abs(signmatrix)
  return(list(adjm=finalmatrix,delay =delaymatrix, cor = cormatrix,sign = signmatrix))
}




