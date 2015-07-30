\name{determinant.lrpd}
\alias{determinant.lrpd}
\title{
Determinant
}
\description{
This function efficiently computes the determinant of an lrpd matrix.
}
\usage{
\method{determinant}{lrpd}(object, logarithm = TRUE, ...)
}
\arguments{
  \item{object}{
an object of class "lrpd".
}
  \item{logarithm}{
logical; If TRUE (default) return the logarithm of the determinant.
}
  \item{\dots}{
not used.
}
}
\value{
\item{scalar}{(logarithm of) the determinant.}
}
\author{
Ye Wang (Eric)

Maintainer: \email{ericwang921198@gmail.com}
}
\seealso{
\code{\link{solve.lrpd}} for inverse.
}
\examples{
library(lrpd)

set.seed(2)
K <- 1000
L <- matrix(rnorm(K*floor(K/10)),K,floor(K/10))
Sl <- matrix(rnorm(floor(K/10)*floor(K/10)),floor(K/10),floor(K/10))
S <- Sl\%*\%t(Sl)+diag(rnorm(floor(K/10))^2)
N <- rnorm(K)^2
R <- L\%*\%S\%*\%t(L) + diag(N)

mat <- lrpd(N,L,S)

system.time(RI1 <- as.numeric(determinant(R)$modulus))
system.time(RI2 <- determinant(mat))

all.equal(RI1,RI2)
}
\keyword{algebra}
