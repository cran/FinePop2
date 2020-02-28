\name{pop_pairwiseFST}
\alias{pop_pairwiseFST}
\title{Genome-wide poppulation-paiwise FST (Nei & Chesser 1983).}

\description{
This function estimates genome-wide poppulation-paiwise FST among subpopulations based on Nei and Chesser's corrected GST (Nei&Chesser 1983) from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}

\usage{
pop_pairwiseFST(popdata)
}

\arguments{
  \item{popdata}{Population data object created by \code{read.GENEPOP} function from a GENEPOP file.}
}

\value{Estimated genome-wide population-pairwise FST.}

\references{
Nei M, Chesser RK (1983) Estimation of fixation indices and gene diversity. \emph{Annals of Human Genetics}, 47, 253-259.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Shuichi Kitada, Hirohisa Kishino}

\seealso{
\code{\link[FinePop2]{read.GENEPOP}},
\code{\link[stats]{as.dist}}, \code{\link[stats]{as.dendrogram}},
\code{\link[stats]{hclust}}, \code{\link[stats]{cmdscale}}
}

\examples{
# Example of GENEPOP file
data(jsmackerel)
jsm.ms.genepop.file <- tempfile()
jsm.popname.file <- tempfile()
cat(jsmackerel$MS.genepop, file=jsm.ms.genepop.file, sep="\n")
cat(jsmackerel$popname, file=jsm.popname.file, sep=" ")

# Data load
# Prepare your GENEPOP file and population name file in the working directory.
# Replace "jsm.ms.genepop.file" and "jsm.popname.file" by your file names.
popdata <- read.GENEPOP(genepop=jsm.ms.genepop.file, popname=jsm.popname.file)

# GST estimation
result.poppairFST <- pop_pairwiseFST(popdata)
poppairFST.d <- as.dist(result.poppairFST)
print(poppairFST.d)

# dendrogram
poppairFST.hc <- hclust(poppairFST.d,method="average")
plot(as.dendrogram(poppairFST.hc), xlab="",ylab="",main="", las=1)

# MDS plot
mds <- cmdscale(poppairFST.d)
plot(mds, type="n", xlab="",ylab="")
text(mds[,1],mds[,2], popdata$pop_names)
}
