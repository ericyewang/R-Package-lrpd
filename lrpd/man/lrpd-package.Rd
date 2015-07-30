\name{lrpd-package}
\alias{lrpd-package}
\docType{package}
\title{
Low-rank Positive Definite Matrices
}
\description{
Efficient matrix inversion and determinant calculation for \emph{low-rank positive definite (lrpd)} matrices. Efficient multivariate Gaussian density calculation and sampling when its covariance matrix is low-rank.
}
\details{
\tabular{ll}{
Package: \tab lrpd\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-07-29\cr
License: \tab GPL-2\cr
}
This pacakge provides efficient \code{solve()} and \code{determinant()} function for class "lrpd". Moreover, the pacakge also provides efficient density calculation, \code{dlrmvnorm()}, and sampling from multivariate Gaussian, \code{rlrmvnorm()}, whose covariance is low-rank.
}
\author{
Ye Wang (Eric)

Maintainer: \email{ericwang921198@gmail.com}
}
\examples{
set.seed(50)
K <- 100
L <- matrix(rnorm(K*floor(K/10)),K,floor(K/10))
Sl <- matrix(rnorm(floor(K/10)*floor(K/10)),floor(K/10),floor(K/10))
S <- Sl\%*\%t(Sl)+diag(rnorm(floor(K/10))^2)
N <- rnorm(K)^2
mat <- lrpd(N,L,S)

# Matrix operation
res1 <- solve(mat)
res2 <- determinant(mat)
# Mutivariate Gaussian
mu <- rnorm(K)
x <- matrix(rnorm(K*10),10,K)
res3 <- dlrmvnorm(x,mu,mat)
res4 <- rlrmvnorm(100,mu,mat)
}
\keyword{array}
\keyword{algebra}
\keyword{multivariate}
