pop_pairwiseFST <-
function(popdata){ #GST Nei & Chesser 1983
  numPop <- popdata$num_pop
  numMarker <- popdata$num_loci
  numInd <- popdata$pop_sizes
  numAllele <- popdata$num_allele
  max.numAllele <- max(numAllele)

  af <- array(0, c(numPop,numMarker,max.numAllele))
  for(cmak in 1:numMarker){
    af[,cmak,1:numAllele[cmak]] <- t(popdata$allele_freq[[cmak]])
  }
  af2 <- af^2

  PopSt <- array(0, c(numPop,numPop))
  dimnames(PopSt) <- list(popdata$pop_names,popdata$pop_names)
  message("Calculating population ", appendLF=FALSE)
  cstep.pop <- ""
  spop <- 2
  for(i in 1:(numPop-1)){
  for(j in (i+1):numPop){
    message(paste0(rep("\b", nchar(cstep.pop)), collapse=""), appendLF=FALSE)
    cstep.pop <- paste0(i, ":", j, " ")
    message(cstep.pop, appendLF=FALSE); flush.console()
    hs <- 1 - colMeans(apply(af2[c(i,j),,], c(1,2), sum))
    ht <- 1 - rowSums(apply(af[c(i,j),,], c(2,3), mean)^2)
    cn <- 1 / mean(1/numInd[c(i,j)])
    Hs <- 2*cn/(2*cn-1) * hs
    Ht <- ht + Hs/(2*cn*spop)
    PopSt[i,j] <- PopSt[j,i] <- 1 - mean(Hs,na.rm=TRUE) / mean(Ht,na.rm=TRUE)
  }}

  message(paste0(rep("\b", nchar(cstep.pop)), collapse=""), appendLF=FALSE)
  message(" done.")

  return(PopSt)
}
