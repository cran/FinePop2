\name{jsmackerel}
\alias{jsmackerel}
\docType{data}
\title{An example dataset of Japanese Spanich mackerel in GENEPOP and frequency format.}

\description{An example of a genetic data for a Japanese Spanish mackerel population (Nakajima et al. 2014). It contains genotypic information of 5 microsatellite markers from 8 subpopulations of 715 individuals. GENEPOP format (Rousset 2008) text files are available. Name list of subpopulations also is attached.}

\usage{data("jsmackerel")}
\format{
 $ MS.genepop: Genotypic information of 5 microsatellites in GENEPOP format text data.\cr
 $ popname: Names of subpopulations.
}

\references{
Nakajima K et al. (2014) Genetic effects of marine stock enhancement: a case study based on the highly piscivorous Japanese Spanish mackerel. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 71, 301-314.\cr\cr
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}

\seealso{
\code{\link[FinePop2]{read.GENEPOP}}
}
\examples{
data(jsmackerel)
jsm.ms.genepop.file <- tempfile()
jsm.popname.file <- tempfile()
cat(jsmackerel$MS.genepop, file=jsm.ms.genepop.file, sep="\n")
cat(jsmackerel$popname, file=jsm.popname.file, sep=" ")

# See two text files in temporary directory.
#  jsm.ms.genepop.file  : GENEPOP format file of microsatellite data
#  jsm.popname.file     : plain text file of subpopulation names
}
\keyword{datasets}
