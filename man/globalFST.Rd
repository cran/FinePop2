\name{globalFST}
\alias{globalFST}
\title{Genome-wide global FST (Weir & Cockerham 1984).}

\description{
This function estimates genom-wide global FST based on Weir and Cockerham's theta (Weir & Cockerham 1984) from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}

\usage{
globalFST(popdata)
}

\arguments{
  \item{popdata}{Population data object created by \code{read.GENEPOP} function from a GENEPOP file.}
}

\value{
\item{fst}{Estimated genome-wide global FST}
\item{se}{Standard error of estimated FST}
}

\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.\cr\cr
Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. \emph{Evolution}, 38, 1358-1370.
}
\author{Reiichiro Nakamichi, Shuichi Kitada, Hirohisa Kishino}


\seealso{
\code{\link[FinePop2]{read.GENEPOP}}
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

# theta estimation
result.globalFST <- globalFST(popdata)
print(result.globalFST)
}
