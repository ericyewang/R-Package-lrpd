\name{solve.lrpd}
\alias{solve.lrpd}
\alias{invlrpd}
\title{
Matrix Inverse
}
\description{
This function efficiently invert an lrpd matrix.
}
\usage{
\method{solve}{lrpd}(a, b, ...)

invlrpd(object,...)
}
\arguments{
  \item{a}{
a square numeric matrix containing the coefficients of the linear system.
}
  \item{b}{
a numeric vector or matrix giving the right-hand side(s) of the linear system. If missing, b is taken to be an identity matrix and solve will return the inverse of a.
}
  \item{object}{
an object of class "lrpd".
}
  \item{\dots}{
not used.
}
}
\value{
\item{matrix}{the solution of the linear system.}
}
\author{
Ye Wang (Eric)

Maintainer: \email{ericwang921198@gmail.com}
}
\seealso{
\code{\link{determinant.lrpd}} for determinant.
}
\examples{
\dontrun{
library(MASS)
library(lrpd)

set.seed(2)
K <- 1000
L <- matrix(rnorm(K*floor(K/10)),K,floor(K/10))
Sl <- matrix(rnorm(floor(K/10)*floor(K/10)),floor(K/10),floor(K/10))
S <- Sl\%*\%t(Sl)+diag(rnorm(floor(K/10))^2)
N <- rnorm(K)^2
R <- L\%*\%S\%*\%t(L) + diag(N)

mat <- lrpd(N,L,S)

# Test matrix inverting
system.time(RI1 <- solve(R))
system.time(RI2 <- chol2inv(chol(R)))
system.time(RI3 <- qr.solve(R))
system.time(RI4 <- solve(mat))

all.equal(RI1, RI2)
all.equal(RI1, RI3)
all.equal(RI1, RI4)
}
}
\keyword{algebra}
