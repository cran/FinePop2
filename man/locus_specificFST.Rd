\name{locus_specificFST}
\alias{locus_specificFST}
\title{Locus-specific global FST (Kitada et al. 2007, 2017).}
\description{
This function estimates locus-specific global FST among subpopulations using empirical Bayes method (Kitada et al. 2007, 2017) from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}

\usage{
locus_specificFST(popdata)
}

\arguments{
  \item{popdata}{Population data object created by \code{read.GENEPOP} function from a GENEPOP file.}
}

\value{Estimated locus-specific global FST.}

\references{
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.\cr\cr
Kitada S, Nakamichi R, Kishino H (2017) The empirical Bayes estimators of fine-scale population structure in high gene flow species. \emph{Mol. Ecol. Resources}, DOI: 10.1111/1755-0998.12663 \cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Shuichi Kitada, Hirohisa Kishino}

\seealso{
\code{\link[FinePop2]{read.GENEPOP}}, \code{\link[FinePop]{read.frequency}}}

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

# FST estimation
locspFST <- locus_specificFST(popdata)
print(locspFST)
}
