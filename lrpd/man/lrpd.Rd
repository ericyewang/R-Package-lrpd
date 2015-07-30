\name{lrpd}
\alias{lrpd}
\alias{print.lrpd}
\alias{summary.lrpd}
\alias{print.summary.lrpd}
\title{
Low-rank Positive Definite Matrices
}
\description{
Run this function to create a "lrpd" object.
}
\usage{
lrpd(N, L, S = diag(rep(1, ncol(L))), ...)

\method{print}{lrpd}(object,...)
\method{summary}{lrpd}(object,...)
}
\arguments{
  \item{N}{
vector containing the diagonal elements.
}
  \item{L}{
matrix with each column being a factor loading vector.
}
  \item{S}{
score matrix.
}
  \item{object}{
object of class "lrpd".
}
  \item{\dots}{
not used.
}
}
\details{
\emph{Low-rank positive definite} matrix M has the following form

\deqn{M=N+LSL^{T},}

where N is a positive diagonal matrix and S is a positive definite matrix.
}
\author{
Ye Wang (Eric)
}
\examples{
set.seed(50)
K <- 100
L <- matrix(rnorm(K*floor(K/10)),K,floor(K/10))
Sl <- matrix(rnorm(floor(K/10)*floor(K/10)),floor(K/10),floor(K/10))
S <- Sl\%*\%t(Sl)+diag(rnorm(floor(K/10))^2)
N <- rnorm(K)^2

mat <- lrpd(N,L,S)
print(mat)
summary(mat)
}
\keyword{array}
