locus_specificFST <-
function(popdata){
nPop <- popdata$num_pop
nLoci <- popdata$num_loci
LocusNames <- popdata$loci_names
PopNames <- popdata$pop_names

nAlleles <- popdata$num_allele
nAllelesMax <- max(nAlleles)
n.lpa <- array(NA, c(nLoci,nPop,nAllelesMax))
for (cl in 1:nLoci){
  obs.an <- t(popdata$allele_count[[cl]])
  check0 <- colSums(obs.an)==0        # remove 0 freq allele
  obs.an <- obs.an[,!check0, drop=FALSE]
  cna <- ncol(obs.an)
  nAlleles[cl] <- cna
  n.lpa[cl,,1:cna] <- obs.an
}
nAllelesMax <- max(nAlleles)
n.lpa <- n.lpa[,,1:nAllelesMax, drop=FALSE]

n.la <- apply(n.lpa, c(1,3), sum, na.rm=TRUE)
n.l <- rowSums(n.la, na.rm=TRUE)
n.lp <- apply(n.lpa, c(1,2), sum, na.rm=TRUE)
n.lp.array <- array(rep(n.lp, nAllelesMax), c(nLoci,nPop,nAllelesMax))

# Excluding one-sided markers (major alelle freq > 99%)
maf_check <- apply(n.la/n.l, 1, max, na.rm=TRUE) > 0.99
if(sum(maf_check)>0){
  message(
    "WARNING: Detected inappropriate markers (major allele frequency > 99%).\n",
    " Locus Names: ", paste(LocusNames[maf_check],collapse=", "))
}

message("Computing global differenciation with ML method... ",appendLF=FALSE);flush.console()
 
p.lpa <- n.lpa / n.lp.array
p.ave.la <- apply(p.lpa, c(1,3), mean, na.rm=TRUE)
p.ave.la.array <- array(rep(p.ave.la, nPop), c(nLoci,nAllelesMax,nPop))
p.ave.la.array <- aperm(p.ave.la.array, c(1,3,2))

sigma2 <- apply((p.lpa - p.ave.la.array)^2, c(1,3), sum, na.rm=TRUE) / (nPop-1)

globalFst.l <- rep(0,nLoci)
nlogL.l <- function(log.theta){
  theta <- exp(log.theta)
  alpha <- theta * cp.ave.la
  alpha.array <- array(rep(alpha,nPop), c(1,nAllelesMax,nPop))
  lgtheta <- lgamma(theta)
  lgalpha <- lgamma(alpha.array)
  clogL <- lgtheta * 1 * nPop - sum(lgamma(cn.lp+theta)) +
           sum(lgamma(cn.lap+alpha.array)-lgalpha, na.rm=TRUE)
  return(-clogL)
}
for(cloc in 1:nLoci){
  cnLoci <- 1
  cp.ave.la <- p.ave.la[cloc,,drop=FALSE]
  cn.lp <- n.lp[cloc,,drop=FALSE]
  theta_init <- mean(cp.ave.la*(1-cp.ave.la)/sigma2[cloc,,drop=FALSE], na.rm=TRUE) - 1
  cn.lap <- aperm(n.lpa[cloc,,,drop=FALSE], c(1,3,2))
  opt_result <- optim(par=log(theta_init), fn=nlogL.l, method="L-BFGS-B",
                      lower=log(theta_init/1000), upper=log(theta_init*1000),
                      hessian=FALSE, control=list(maxit=1000))
  ctheta_est <- exp(opt_result$par)
  globalFst.l[cloc] <- 1/(ctheta_est+1)
}
message("done.")
return(globalFst.l)
}
