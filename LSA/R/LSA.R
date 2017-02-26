#Part I - LSA Code
#i: represents varying on unit i
#j: represents varying on unit j
#w: represents varying on timeslot w
#d: represents varying on delay
#Order of subscripts indicates order in array, row, column... etc..

#POPiw = population for unit i at timeslot k
#RNiw = rank normal score for unit i at timeslot k

#LSAij_d = Local similarity score for unit i, unit j, delay d
#CORNij_d = correlation using rank normal for unit i, unit j, delay d
#CORij_d = correlation for unit i, unit j, delay d

#CORNSij_d = Correlation (using rank normal) significance for unit i, unit j, delay d
#CORSij_d = Correlation significance for unit i, unit j, delay d


###Note: I could split this up into 2 timeslots to track the exact pair of times that produced a particular
###Combination, but as of yet, I can't use this, better to make the array more manageable and just track 1 slot
###All I need right now is the series for a given delay, not all the exact times

#Tw = vector of absolute times associated with timeslot k

#returns LSA score
LSA <- function(x) {
  yP = 0
  yN = 0
  yMaxP = 0
  yMaxN = 0
  tmp = x
  tmp[is.na(tmp)]=0

  for(i in tmp) {
    yP = max(yP+i,0)
    yN = min(yN+i,0)
    yMaxP = max(yMaxP,yP)
    yMaxN = min(yMaxN,yN)
  }
  if(yMaxP>abs(yMaxN)) {
    return(yMaxP)
  }
  else {
    return(yMaxN)
  }

}


normalTransform <- function(POPiw) t(apply(POPiw,1,function(x) qnorm(rank(x,na.last="keep")/(1+sum(!is.na(x))))))

#widens Tk to account for gaps
expandTw <- function(Tw,res=1) seq(min(Tw),max(Tw),by=res)

#widens Matrix according to Tk
expandMatrixUsingTw <- function(timeMatrix,Tw,res=1) {
  TwE <- expandTw(Tw,res=res)
  timeMatrixNew <- matrix(NA,nrow=nrow(timeMatrix),ncol=length(TwE))
  timeMatrixNew[,is.element(TwE,Tw)] = timeMatrix
  return(timeMatrixNew)
}


#res should always divide maxdelay and the minimum time measurement
compute_LSA_COR <- function(POPiw,maxdelay = 1,res=1,Tw=c(),expand=!is.null(Tw)) {
  POPiw_ = POPiw
  Tw_ = Tw

  #expand will expand the matrix given the time.vector to ensure resolution is standard
  if(expand) {
    POPiw_=expandMatrixUsingTw(POPiw,Tw,res)
    Tw_ = expandTw(Tw)
  }

  #number of OTUS
  n = dim(POPiw_)[[1]]
  #number of timeslots
  m = dim(POPiw_)[[2]]

  RNiw = normalTransform(POPiw_)

  #Store the delay range
  delrange = seq(0,maxdelay,by=res)
  ndelay = length(delrange)

  LSAij_d = array(NA,dim=c(n,n,ndelay))
  dimnames(LSAij_d)[[1]] = paste("OTU",1:n,sep="")
  dimnames(LSAij_d)[[2]] = paste("OTU",1:n,sep="")
  dimnames(LSAij_d)[[3]] = paste("DELAY",delrange, sep="")

  CORNij_d = array(NA,dim=c(n,n,ndelay))
  dimnames(CORNij_d)[[1]] = paste("OTU",1:n,sep="")
  dimnames(CORNij_d)[[2]] = paste("OTU",1:n,sep="")
  dimnames(CORNij_d)[[3]] = paste("DELAY",delrange, sep="")

  CORij_d = array(NA,dim=c(n,n,ndelay))
  dimnames(CORij_d)[[1]] = paste("OTU",1:n,sep="")
  dimnames(CORij_d)[[2]] = paste("OTU",1:n,sep="")
  dimnames(CORij_d)[[3]] = paste("DELAY",delrange, sep="")

  CORNSij_d = array(NA,dim=c(n,n,ndelay))
  dimnames(CORNSij_d)[[1]] = paste("OTU",1:n,sep="")
  dimnames(CORNSij_d)[[2]] = paste("OTU",1:n,sep="")
  dimnames(CORNSij_d)[[3]] = paste("DELAY",delrange, sep="")

  CORSij_d = array(NA,dim=c(n,n,ndelay))
  dimnames(CORSij_d)[[1]] = paste("OTU",1:n,sep="")
  dimnames(CORSij_d)[[2]] = paste("OTU",1:n,sep="")
  dimnames(CORSij_d)[[3]] = paste("DELAY",delrange, sep="")

  for(k in 1:ndelay) {

    #We will only be taking part of the time series range
    #to account for this, we need to remove data points when considering a delay
    #in this case, positive delay means that the second is offset from the beginning,
    #and the first is offset from the end
    d = delrange[k]
    wrange = 1:(m-d/res)
    vrange = (1+d/res):m

    #compute LSA and COR on every pair of rows and store by delay

    LSAij = t(apply(RNiw[,wrange],1,function(x) apply(RNiw[,vrange],1,function(y) LSA(x*y))))


    CORNij = t(apply(RNiw[,wrange],1,function(x) apply(RNiw[,vrange],1,function(y) {
      not_na = !is.na(x) & !is.na(y)
      if(sum(not_na)>10) return(cor.test(x[not_na], y[not_na])$estimate)
      else return(NA)
    })))


    CORij = t(apply(POPiw_[,wrange],1,function(x) apply(POPiw_[,vrange],1,function(y) {
      not_na = !is.na(x) & !is.na(y)
      if(sum(not_na)>10) return(cor.test(x[not_na], y[not_na])$estimate)
      else return(NA)

    })))

    CORNSij = t(apply(RNiw[,wrange],1,function(x) apply(RNiw[,vrange],1,function(y) {
      not_na = !is.na(x) & !is.na(y)
      if(sum(not_na)>10) return(cor.test(x[not_na], y[not_na])$p.value)
      else return(NA)

    })))

    CORSij = t(apply(POPiw_[,wrange],1,function(x) apply(POPiw_[,vrange],1,function(y) {
      not_na = !is.na(x) & !is.na(y)
      if(sum(not_na)>10) return(cor.test(x[not_na], y[not_na])$p.value)
      else return(NA)
    })))

    dim(LSAij) = c(n,n)
    dim(CORNij) = c(n,n)
    dim(CORij) = c(n,n)
    dim(CORij) = c(n,n)
    dim(CORSij) = c(n,n)

    LSAij_d[,,k] = LSAij
    CORNij_d[,,k] = CORNij
    CORij_d[,,k] = CORij
    CORNSij_d[,,k] = CORNSij
    CORSij_d[,,k] = CORSij

  }





  return(list(lsa_d = LSAij_d,
              corn_d = CORNij_d,
              cor_d = CORij_d,
              corns_d = CORNSij_d,
              cors_d = CORSij_d
  ))
}

ConfidentPairs <- function(LSA_COR_obj,lsasighash_file="LSA98hashtable.csv",origdata = m000.matrix[OTUS,],cinterv=0.95) {
  lsadata = LSA_COR_obj$lsa_d
  corsdata = LSA_COR_obj$cors_d
  cornsdata = LSA_COR_obj$corns_d
  hashtable = read.csv(lsasighash_file)
  datazeros = apply(origdata,1,function(x) sum(x==0,na.rm=TRUE))

  dim1 = dim(lsadata)


  #significance level is determined by the number of ties in the original data
  #we expect these ties to take place at 0, effecting the significance value of the score
  #hence a full set of hashed permutations of the data is stored prior to the analysis
  #this strongly cuts down on the number of duplicate permutation runs
  #one way to change this would be to ensure that all vectors have rank score conserved
  #but this may bias certain results

  lsasiglevels = matrix(Inf,nrow=dim1,ncol=dim1)
  for(i in 1:(dim1[[1]]-1)) for(j in (i+1):dim1[[2]]) lsasiglevels[i,j] = hashtable[datazeros[i],datazeros[j]]

  #compute COR sig for i j as the minimum over all delays
  COR_strong = apply(corsdata,c(1,2),function(x) min(x,na.rm=TRUE))

  #will contain numeric(0)
  COR_delay = apply(corsdata,c(1,2), function(x) which.max(x)-1)

  #compute CORN sig for i j as the minimum over all delays
  CORN_strong = apply(cornsdata,c(1,2),function(x) min(x,na.rm=TRUE))

  #will contain numeric(0)
  CORN_delay = apply(cornsdata,c(1,2), function(x) which.max(x)-1)

  #compute LSA for i j as the maximum over all delays
  LSA_strong = apply(lsadata,c(1,2),function(x) max(abs(x),na.rm=TRUE))

  #need to fix this so that it gives the proper delays, but currently, this should be alright
  LSA_delay = apply(lsadata,c(1,2),function(x) which.max(abs(x))-1)


  sigLSA = LSA_strong>lsasiglevels
  sigCOR = COR_strong<(1-cinterv)
  sigCORN = CORN_strong<(1-cinterv)

  sigLSA[is.na(sigLSA)]=FALSE
  sigCOR[is.na(sigCOR)]=FALSE
  sigCORN[is.na(sigCORN)]=FALSE

  #don't want the self-score or duplicate scores
  selfscore = matrix(FALSE,dim1[[1]],dim1[[1]])
  for(i in 1:dim1[[1]]) for(j in i:dim1[[1]]) selfscore[j,i]=TRUE


  lsapairs = which(sigLSA & !selfscore,arr.ind=TRUE)
  corpairs = which(sigCOR & !selfscore,arr.ind=TRUE)
  cornpairs = which(sigCORN & !selfscore,arr.ind=TRUE)
  allpairs = which((sigLSA | sigCOR | sigCORN) & !selfscore,arr.ind=TRUE)
  sharedpairs = which(sigLSA & sigCOR & sigCORN & !selfscore,arr.ind=TRUE)

  return(list(lsa_pairs=lsapairs,cor_pairs=corpairs,corn_pairs=cornpairs,all_pairs=allpairs,shared_pairs=sharedpairs))

}

pairlistToCSV <- function(pairlist,filename="LSACORpairs.csv") {
  pairs = pairlist$all_pairs
  totpairs = dim(pairs)[[1]]
  csvoutput = matrix(ncol=3,nrow=totpairs)
  csvoutput[,1] = pairs[,1]
  csvoutput[,2] = "pu"
  csvoutput[,3] = pairs[,2]
  colnames(csvoutput) = c("Source","Interaction","Target")
  write.csv(csvoutput,file = filename,quote=FALSE)
  return(filename)
}




