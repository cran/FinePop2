\name{pop_specificFST}
\alias{pop_specificFST}
\title{Genome-wide population-specific FST (Weir & Goudet 2017).}

\description{
This function estimates genome-wide poppulation-specific FST based on Weir and Goudet's Method (Weir&Goudet 2017) from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}

\usage{
pop_specificFST(popdata, cov = FALSE)
}

\arguments{
  \item{popdata}{Population data object created by \code{read.GENEPOP} function from a GENEPOP file.}
  \item{cov}{A Logical argument indicating whether variance-covariance matrix of estimated FST should be calculated.}
}

\value{
\item{fst}{Estimated genome-wide population-specific FST and standard error.}
\item{cov}{Variance-covariance matrix of estimated FST}
}

\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.\cr\cr
Weir BS, Goudet J (2017) A Unified Characterization of Population Structure and Relatedness. \emph{Genetics}, 206, 2085-2103.
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

# FST estimation
result.popspFST <- pop_specificFST(popdata)
print(result.popspFST$fst)
}
