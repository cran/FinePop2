read.GENEPOP <-
function(genepop, popname=NULL){
# read genepop file
all_lines <- scan(genepop, what=character(), quiet=TRUE, sep="\n", blank.lines.skip=FALSE)

# title 
gp_title <- all_lines[1]

# locus pop sample count
cline <- gsub(" ", "", all_lines)
cline <- gsub("\t", "", cline)
poploc <- which(toupper(cline)=="POP")
MarkerList <- cline[2:(poploc[1]-1)]
MarkerList <- gsub(",", "\b", MarkerList)
MarkerList <- unlist(strsplit(MarkerList, "\n"))
numMarker <- length(MarkerList)
numPop <- length(poploc)
popstart <- poploc + 1
popend <- c(poploc[-1]-1, length(all_lines))
numInd <- popend - popstart + 1
numIndAll <- sum(numInd)
PopID <- rep(1:numPop, numInd)
rm(cline);gc()

if(is.null(popname)){
  pop.names <- paste0("pop",1:numPop)
}else{
  pop.names <- scan(popname, what=character(), quiet=TRUE, blank.lines.skip=TRUE)
}

# marker genotype
gtdata <- all_lines[-poploc]
gtdata <- gtdata[-(1:(poploc[1]-1))]
gtdata <- unlist(strsplit(gtdata, ","))
IndID <- gtdata[c(TRUE,FALSE)]
IndID <- gsub(" ", "", IndID)
IndID <- gsub("\t", "", IndID)
gtdata <- gsub(" ", "\t", gtdata[c(FALSE,TRUE)])
gtdata <- matrix(unlist(strsplit(gtdata, "\t")), nrow=numIndAll, byrow=TRUE)
gtdata <- gtdata[,-which(colSums(gtdata=="")!=0)]
gp_digit <- as.integer(nchar(as.character(gtdata[1,1]))/2)
gp_na <- paste(rep("0", gp_digit), collapse="")
rm(all_lines);gc()

# gpdata
htdata1 <- substr(gtdata,1,gp_digit)
htdata2 <- substr(gtdata,gp_digit+1,gp_digit*2)
htdata1 <- gsub(gp_na, NA, htdata1)
htdata2 <- gsub(gp_na, NA, htdata2)

Haplotype <- list()
Ho <- rep(0, numPop)
for(cpop in 1:numPop){
  cpopind <- PopID==cpop
  cnpop <- numInd[cpop]
  chaplo <- array("",c(cnpop,numMarker,2))
  h1 <- htdata1[cpopind,,drop=FALSE]
  h2 <- htdata2[cpopind,,drop=FALSE]
  chaplo[,,1] <- h1
  chaplo[,,2] <- h2
  Ho[cpop] <- mean(colMeans(h1!=h2,na.rm=TRUE),na.rm=TRUE)
  Haplotype[[cpop]] <- chaplo
}

call <- !is.na(htdata1)
CallRate.loci <- colMeans(call)
CallRate.ind <- list()
call.ind <- rowMeans(call)
for(cp in 1:numPop){CallRate.ind[[cp]] <-call.ind[PopID==cp]}

htdata <- rbind(htdata1, htdata2)
rm(htdata1, htdata2, gtdata);gc()

AlleleCount <- list()
AlleleFreq <- list()
IndObs <- list()
numAlleles <- rep(0,numMarker)
AlleleList <- list()
for(cm in 1:numMarker){
  cgt <- table(htdata[,cm], c(PopID, PopID), useNA="no")
  colnames(cgt) <- NULL
  numAlleles[cm] <- nrow(cgt)
  AlleleList[[cm]] <- rownames(cgt)
  cgtnum <- colSums(cgt)
  AlleleCount[[cm]] <- cgt
  numcall <- as.integer(cgtnum/2)
  IndObs[[cm]] <- numcall
  AlleleFreq[[cm]] <- t(t(cgt) / cgtnum) 
}
He <- sapply(AlleleFreq, function(x){1 - colSums(x^2)})
He <- rowMeans(He,na.rm=TRUE)

IndNames <- list()
for(cpop in 1:numPop){IndNames[[cpop]] <- IndID[PopID==cpop]}

return(list(genotype=Haplotype,
            allele_count=AlleleCount,
            allele_freq=AlleleFreq,
            ind_count=IndObs,
            num_pop=numPop,
            pop_sizes=numInd,
            pop_names=pop.names,
            ind_names=IndNames,
            num_loci=numMarker,
            loci_names=MarkerList,
            num_allele=numAlleles,
            allele_list=AlleleList,
            call_rate_loci=CallRate.loci,
            call_rate_ind=CallRate.ind,
            He=He,
            Ho=Ho
       ))
}
