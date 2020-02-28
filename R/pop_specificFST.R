pop_specificFST <-
function(popdata, cov=FALSE){ #pop-sp FST Weir & Goudet 2017
  numpop <- popdata$num_pop
  numloci <- popdata$num_loci
  ni <- popdata$pop_sizes

  # FST
  fst.wg.num <- array(0, c(numloci,numpop))
  fst.wg.den <- rep(0, numloci)
  message("Calculating population ", appendLF=FALSE)
  cstep.pop <- ""
  for(cl in 1:numloci){
    message(paste0(rep("\b", nchar(cstep.pop)), collapse=""), appendLF=FALSE)
    cstep.pop <- paste0(cl, "/", numpop)
    caf.pop <- popdata$allele_freq[[cl]]
    cni <- popdata$ind_count[[cl]]
    cMwi <- 2*cni/(2*cni-1) * colSums(caf.pop^2,na.rm=TRUE) - 1/(2*cni-1)
    cMBii <- 0
    for(cp1 in 1:(numpop-1)){
    for(cp2 in (cp1+1):numpop){
      cMBii <- cMBii + sum(caf.pop[,cp1] * caf.pop[,cp2], na.rm=TRUE)
    }}
    cMBii <- cMBii * 2
    cMB <- 1/(numpop*(numpop-1)) * cMBii
    fst.wg.num[cl,] <- cMwi - cMB
    fst.wg.den[cl] <- 1 - cMB
  }
  popfst.wg <- colMeans(fst.wg.num) / mean(fst.wg.den)
  names(popfst.wg) <- popdata$pop_names

  # Variance/Covariance of FST
  xbar <- mean(fst.wg.den)
  ybar <- colMeans(fst.wg.num)
  V_popfst.wg <- (ybar^2)/(xbar^4)*var(fst.wg.den)/numloci +
                   1/(xbar^2) * apply(fst.wg.num,2,var)/numloci -
                   (2*ybar)/(xbar^3) * as.numeric(cov(fst.wg.num,fst.wg.den))/numloci

  COV_popfst.wg <- NULL
  if(cov){
    COV_popfst.wg <- array(NA, c(numpop,numpop))
    for(cp1 in 1:(numpop-1)){
      ybar1 <- ybar[cp1]
    for(cp2 in (cp1+1):numpop){
      ybar2 <- ybar[cp2]
      COV_popfst.wg[cp1,cp2] <- COV_popfst.wg[cp2,cp1] <-
        ybar1*ybar2/(xbar^4)*var(fst.wg.den)/numloci -
        ybar1/(xbar^3) * cov(fst.wg.num[,cp1],fst.wg.den)/numloci -
        ybar2/(xbar^3) * cov(fst.wg.num[,cp2],fst.wg.den)/numloci
    }}
    diag(COV_popfst.wg) <- V_popfst.wg
  }

  message(paste0(rep("\b", nchar(cstep.pop)), collapse=""), appendLF=FALSE)
  message(" done.")

  return(list(
    fst=data.frame(FST=popfst.wg, SE=sqrt(V_popfst.wg)),
    cov=COV_popfst.wg
  ))
}
