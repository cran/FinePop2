\name{herring}
\alias{herring}
\docType{data}
\title{An example dataset of Atlantic herring.}
\description{
An example of a genetic data for Atlantic herring population (Limborg et al. 2012). It contains genotypic information of 281 SNPs from 18 subpopulations of 607 individuals. GENEPOP format (Rousset 2008) text file is available. Subpopulation names, environmental factors (longitude, latitude, temperature and salinity) at each subpopulation are attached.}
\usage{data("herring")}
\format{
 $ genepop       : Genotypic information of 281 SNPs in GENEPOP format text data.\cr
 $ popname       : Names of subpopulations.\cr
 $ environment   : Table of temperature and salinity at each subpopulation.\cr
}
\references{
Limborg MT, Helyar SJ, de Bruyn M et al. (2012) Environmental selection on transcriptome-derived SNPs in a high gene flow marine fish, the Atlantic herring (Clupea harengus). \emph{Molecular Ecology}, 21, 3686-3703. \cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\seealso{
\code{\link[FinePop2]{read.GENEPOP}}, \code{\link[FinePop2]{GLS}}
}
\examples{
data(herring)
ah.genepop.file <- tempfile()
ah.popname.file <- tempfile()
cat(herring$genepop, file=ah.genepop.file, sep="\n")
cat(herring$popname, file=ah.popname.file, sep=" ")

# See two text files in temporary directory.
#  ah.genepop.file  : GENEPOP format file of 281SNPs in 18 subpopulations
#  ah.popname.file  : plain text file of subpopulation names

print(herring$environment)
}
\keyword{datasets}
