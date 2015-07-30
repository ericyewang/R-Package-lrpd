\name{lrmvnorm}
\alias{lrmvnorm}
\alias{dlrmvnorm}
\alias{rlrmvnorm}
\title{
Multivariate Gaussian with Low-rank Covariance
}
\description{
Efficient density calculation and sampling.
}
\usage{
dlrmvnorm(x, mu, Sigma, logarithm = TRUE, ...)
rlrmvnorm(n, mu, Sigma, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
}
  \item{n}{
number of observations.
}
  \item{mu}{
mean vector.
}
  \item{Sigma}{
covariance matrix.
}
  \item{logarithm}{
logical; if TRUE (default), the logarithm of the densities are given.
}
  \item{\dots}{
not used.
}
}
\author{
Ye Wang (Eric)

Maintainer: \email{ericwang921198@gmail.com}
}
\examples{
\dontrun{
library(mvtnorm)
library(lrpd)

set.seed(2)
K <- 1000
L <- matrix(rnorm(K*floor(K/10)),K,floor(K/10))
Sl <- matrix(rnorm(floor(K/10)*floor(K/10)),floor(K/10),floor(K/10))
S <- Sl\%*\%t(Sl)+diag(rnorm(floor(K/10))^2)
N <- rnorm(K)^2
R <- L\%*\%S\%*\%t(L) + diag(N)
# Gaussian mean vectors and synthethic data
mu <- rnorm(K)
x <- matrix(rnorm(K*10),10,K)

mat <- lrpd(N,L,S)

# logarithm of densities
system.time(RI1 <- dmvnorm(x,mu,R,log=TRUE))
system.time(RI2 <- dlrmvnorm(x,mu,mat))

all.equal(RI1,RI2)

# random samples
system.time(RI1 <- rmvnorm(100,mu,R))
system.time(RI2 <- rlrmvnorm(100,mu,mat))
}
}
\keyword{multivariate}
