\name{mult_diag}
\alias{mult_diag}
\title{
Efficient Matrix Multiplication with a Diagonal Matrix
}
\description{
This functoin efficiently computes the matrix multiplication between an arbitrary matrix and a diagonal matrix.
}
\usage{
mult_diag(a, b, ...)
}
\arguments{
  \item{a}{
matrix to be multiplied with a diagonal matrix.
}
  \item{b}{
vector containing the diagonal elements of the diagonal matrix; if length(b)==ncol(a), the function returns a%*%diag(b), if length(b)==nrow(a), the function returns diag(b)%*%a, otherwise the function will generate an error message.
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
n1 <- 1000; n2 <- 100
A <- matrix(rnorm(n1*n2),n1,n2)
d <- rnorm(n1)
system.time(R1 <- diag(d)%*%A)
system.time(R2 <- mult_diag(A,d))

all.equal(R1, R2)
}
}
\keyword{algebra}
