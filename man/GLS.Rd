\name{GLS}
\alias{GLS}
\title{Generalized least squres for regression analysis considering auto correlation.}

\description{
This function provides a multiple regression analysis considering auto correlation of response variable using generalized least squres method (Aitken 1934). It supports \code{\link[stats]{lm}} like format of \code{model}.  A typical model has the form \code{response ~ terms}. Terms specification supports only \code{first + second} form. Cross term specification of \code{first * second} form is not supported.
}

\usage{
GLS(model, data, omega = NULL)
}

\arguments{
  \item{model}{Symbolic description of the model to be fitted.}
  \item{data}{Data frame containing variables in the \code{model}.}
  \item{omega}{A numeric matrix of auto correlation of responce variable.}
}

\value{
\item{coefficients}{Estimated coefficient, standard error, Z value and p value of each factor.}
\item{variance}{Variance-covariance matrix of estimated coefficients.}
\item{logL}{Log likelihood of fitted model.}
}
\references{
Aitken AC (1934) On Least-squares and Linear Combinations of Observations. \emph{Proceedings of the Royal Society of Edinburgh}, 55, 42-48.}
\author{Reiichiro Nakamichi, Shuichi Kitada, Hirohisa Kishino}

\seealso{
\code{\link[stats]{lm}}
}
\examples{
# Example data of Atlantic herring
data(herring)
ah.genepop.file <- tempfile()
ah.popname.file <- tempfile()
cat(herring$genepop, file=ah.genepop.file, sep="\n")
cat(herring$popname, file=ah.popname.file, sep=" ")

# Data load
popdata <- read.GENEPOP(ah.genepop.file, ah.popname.file)

# Pop-specific FST and correlation among populations
fst.popsp <- pop_specificFST(popdata, cov=TRUE)
cov.fst.popsp <- fst.popsp$cov
sd.fst.popsp <- sqrt(diag(cov.fst.popsp))
cov2.fst.popsp <- apply(cov.fst.popsp, 2, function(x){x / sd.fst.popsp})
cor.fst.popsp <- apply(cov2.fst.popsp, 1, function(x){x / sd.fst.popsp})

# Pop-pairwise FST and population structure
fst.poppair <- pop_pairwiseFST(popdata)
fst.md <- cmdscale(fst.poppair)

# GLS analysis of FST and environmental factors
test.data <- data.frame(fst=fst.md[,1], herring$environment)
GLS(fst~., scale(test.data), omega=cor.fst.popsp)
}
