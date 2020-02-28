\name{read.GENEPOP}
\alias{read.GENEPOP}
\title{Create a genotype data object of populations from a GENEPOP format file.}

\description{
This function reads a GENEPOP format file (Rousset 2008) and parse it into an R data object. This data object provides a summary of genotype/haplotype of each sample, allele frequency in each population, and marker status. This data object is used in downstream analysis of this package.
}

\usage{
read.GENEPOP(genepop, popname = NULL)
}

\arguments{
  \item{genepop}{A character value specifying the name of the GENEPOP file to be analyzed.}
  \item{popname}{A character value specifying the name of the plain text file containing the names of subpopulations to be analyzed. This text file must not contain other than subpopulation names. The names must be separated by spaces, tabs or line breaks. If this argument is omitted, serial numbers will be assigned as subpopulation names.}
}

\value{
  \item{num_pop}{Number of subpopulations.}
  \item{pop_sizes}{Number of samples in each subpopulation.}
  \item{pop_names}{Names of subpopulations.}
  \item{ind_names}{Names of samples in each subpopulation.}
  \item{num_loci}{Number of loci.}
  \item{loci_names}{Names of loci.}
  \item{num_allele}{Number of alleles at each locus.}
  \item{allele_list}{A list of alleles at each locus.}
  \item{ind_count}{Observed count of genotyped samples in each subpopulation at each locus.}
  \item{allele_count}{Observed count of genotyped alleles in each subpopulation at each locus.}
  \item{allele_freq}{Observed allele frequencies in each subpopulation at each locus.}
  \item{genotype}{Genotypes of each sample at each locus in haploid designation.}
  \item{call_rate_loci}{Call rate of each locus (rate of genotyped samples at each locus).}
  \item{call_rate_ind}{Call rate of each sample (rate of genotyped markers for each sample).}
  \item{He}{Expected heterozigosity in each subpopulation.}
  \item{Ho}{Observed heterozigosity in each subpopulation.}
}

\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi}

\examples{
# Example of GENEPOP file
data(jsmackerel)
jsm.ms.genepop.file <- tempfile()
jsm.popname.file <- tempfile()
cat(jsmackerel$MS.genepop, file=jsm.ms.genepop.file, sep="\n")
cat(jsmackerel$popname, file=jsm.popname.file, sep=" ")

# Read GENEPOP file with subpopulation names.
# Prepare your GENEPOP file and population name file in the working directory.
# Replace "jsm.ms.genepop.file" and "jsm.popname.file" by your file names.
popdata <- read.GENEPOP(genepop=jsm.ms.genepop.file, popname=jsm.popname.file)

# Read GENEPOP file without subpopulation names.
popdata.noname <- read.GENEPOP(genepop=jsm.ms.genepop.file)
}
