globalFST <-
function(popdata){
npop <- popdata$num_pop
nlocus <- popdata$num_loci
nsamples <- popdata$pop_sizes

## global FST ##
message("Calculating global FST ... ", appendLF=FALSE);flush.console()
sum_a <- sum_abc <- 0
sum_a.l <- sum_abc.l <- rep(0, nlocus)
for(cloc in 1:nlocus){
  nalleles <- popdata$num_allele[cloc]
  allele_names <- rownames(popdata$allele_freq[[cloc]])
  n.obsamples <- popdata$ind_count[[cloc]]

  genotypes <- array(NA, c(npop,max(nsamples),2))
  for(cpop in 1:npop){
    genotypes[cpop,1:nsamples[cpop],1] <- popdata$genotype[[cpop]][,cloc,1]
    genotypes[cpop,1:nsamples[cpop],2] <- popdata$genotype[[cpop]][,cloc,2]
  }

  csum_a.l <- csum_abc.l <- 0
  for(callele in 1:nalleles){
    n_bar <- mean(n.obsamples)
    n_c <- (n_bar * npop - sum(n.obsamples^2)/(n_bar*npop)) / (npop-1)

    freqA <- popdata$allele_freq[[cloc]][callele,]
    p_bar <- sum(freqA * n.obsamples) / (npop*n_bar)
    s2 <- sum(n.obsamples * (freqA - p_bar)^2) / ((npop-1)*n_bar)

    callelename <- allele_names[callele]
    heteros <- apply(genotypes, c(1,2), function(x){
      ((x[1]==callelename)&(x[2]!=callelename)) |
      ((x[1]!=callelename)&(x[2]==callelename))
    })

    freqAH <- rowSums(heteros,na.rm=TRUE) / n.obsamples
    h_bar <- sum(freqAH * n.obsamples) / (npop*n_bar)
    WCa <- n_bar/n_c * (s2 - 1/(n_bar-1) * (p_bar*(1-p_bar) - (npop-1)/npop * s2 - h_bar/4 ))
    WCb <- n_bar/(n_bar-1) * (p_bar*(1-p_bar) - (npop-1)/npop * s2 - (2*n_bar-1)/(4*n_bar) * h_bar)
    WCc <- h_bar/2
    if(is.finite(WCa)){
       sum_a <- sum_a + WCa;
       sum_abc <- sum_abc + WCa + WCb + WCc
       csum_a.l <- csum_a.l + WCa;
       csum_abc.l <- csum_abc.l + WCa + WCb + WCc
    }
  }#allele
  sum_a.l[cloc] <- csum_a.l
  sum_abc.l[cloc] <- csum_abc.l 
}#loc
theta.g <- sum_a / sum_abc

## Variance of FST ##
xbar <- mean(sum_abc.l)
ybar <- mean(sum_a.l)
V <- (ybar^2)/(xbar^4)*var(sum_abc.l)/nlocus +
     1/(xbar^2) * var(sum_a.l)/nlocus -
     (2*ybar)/(xbar^3) * cov(sum_a.l,sum_abc.l)/nlocus

message("done.")
return(list(fst=theta.g, se=sqrt(V)))
}
