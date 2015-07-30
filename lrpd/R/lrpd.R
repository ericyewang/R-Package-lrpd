lrpd <- function(N,L,S=diag(rep(1,ncol(L))),...){
  # Construct a low-rank positive definite matrix object
  # Input:
  #   N: Diagonal noises (vectors);
  #   L: Factor loadings;
  #   S: Loading scores (matrix).
  # Output:
  #   res: an object of class "lrpd"
  
  S <- as.matrix(S)
  L <- as.matrix(L)
  N <- as.numeric(N)
  # Error handling
  if (length(N)!=nrow(L)){
    stop("Dimension mismatch: N and L!")
  }
  if (ncol(S)!=ncol(L)){
    stop("Dimension mismatch: S and L!")
  }
  
  res <- list(N = N, 
              L = L, 
              S = S)
  class(res) <- "lrpd"
  res
}

print.lrpd <- function(object,...){
  mat <- diag(object$N) + object$L%*%object$S%*%t(object$L)
  print(as.matrix(mat))
}

summary.lrpd <- function(object,...){
  d <- length(object$N)
  p <- ncol(object$L)
  if (sum(abs(object$N-mean(object$N)))==0){
    homo <- "Yes"
  }else{
    homo <- "No"
  }
  res <- list(d = d,
              p = p,
              homo = homo)
  class(res) <- "summary.lrpd"
  res
}

print.summary.lrpd <- function(object,...){
  cat("Ambient dimension:\n")
  print(object$d)
  cat("Intrinsic dimension:\n")
  print(object$p)
  cat("Homogeneous noise:\n")
  print(object$homo)
}

invlrpd <- function(object,...){
  # Inverse of a lrpd object
  # Input:
  #   object: an lrpd object.
  # Output:
  #   res: the inverse matrix.
  p <- ncol(object$L)
  S <- diag(rep(1,p))+object$S%*%t(object$L)%*%mult_diag(object$L,1/object$N)
  
  LN <- mult_diag(object$L,1/object$N)
  diag(1/object$N)-LN%*%solve(S)%*%object$S%*%t(LN)
}

solve.lrpd <- function(a,b,...){
  if (missing(b)){
    invlrpd(a)
  }else{
    b.tmp <- as.matrix(b)
    res <- invlrpd(a)%*%b.tmp
  }
}

determinant.lrpd <- function(object,logarithm = TRUE,...){
  p <- ncol(object$L)
  S <- diag(rep(1,p)) + object$S%*%t(object$L)%*%mult_diag(object$L,1/object$N)
  Sdet <- determinant(S)$modulus
  res <- as.numeric(sum(log(object$N)) + Sdet)
  if (logarithm){
    res
  }else{
    exp(res)
  }
}

mult_diag <- function(a,b,...){
  # Mutiply a matrix with a diagonal matrix
  # Input:
  #   a: a full matrix;
  #   b: diagonal elements of the diagonal matrix
  # Output:
  #   res: resulted matrix.
  a <- as.matrix(a)
  b <- as.numeric(b)
  dr <- nrow(a)
  dc <- ncol(a)
  if (length(b)==1){
    b <- rep(b,dc)
    apply(as.matrix(1:dc),1,FUN=function(x){
      return(a[,x]*b[x])
    })
  }else{
    if (length(b)==dc){
      apply(as.matrix(1:dc),1,FUN=function(x){
        return(a[,x]*b[x])
      })
    }else if (length(b)==dr){
      t(apply(as.matrix(1:dr),1,FUN=function(x){
        return(a[x,]*b[x])
      }))
    }else{
      stop("Dimension mismatch!")
    }
  }
}

dlrmvnorm <- function(x,mu,Sigma,logarithm = TRUE,...){
  # density of low-rank multivariate Gaussian
  # Input:
  #   x: data vector (or matrix with each row being a data vector);
  #   mu: mu vector;
  #   Sigma: lrpd object.
  # Output:
  #   res: density.
  
  d <- length(Sigma$N)
  mu <- as.numeric(mu)
  # Error handling
  if (length(mu)!=d){
    stop("Dimension mismatch: mu and Sigma!")
  }
  
  if (is.vector(x)){
    tmp <- t(x-mu)
    n <- 1
  }else{
    n <- nrow(x)
    # error handling
    if (ncol(x)!=d){
      stop("Dimenison mismatch: x and Sigma!")
    }
    tmp <- t(apply(x,1,FUN=function(x){
      x-mu
    }))
  }
  ld <- -log(2*pi)*d/2 - determinant(Sigma)/2
  p <- ncol(Sigma$L)
  S <- diag(rep(1,p))+Sigma$S%*%t(Sigma$L)%*%mult_diag(Sigma$L,1/Sigma$N)
  Sinv <- solve(S)
  ldtmp <- as.numeric(apply(t(1:n),2,FUN=function(x){
    Ltmp <- tmp[x,]%*%mult_diag(Sigma$L,1/Sigma$N)
    ld - sum(tmp[x,]^2/Sigma$N)/2 + Ltmp%*%Sinv%*%Sigma$S%*%t(Ltmp)/2
  }))
  if (logarithm){
    ldtmp
  }else{
    exp(ldtmp)
  }
}

rlrmvnorm <- function(n,mu,Sigma,...){
  # density of low-rank multivariate Gaussian
  # Input:
  #   n: sample size;
  #   mu: mu vector;
  #   Sigma: lrpd object.
  # Output:
  #   res: density.
  
  d <- length(Sigma$N)
  p <- ncol(Sigma$L)
  mu <- as.numeric(mu)
  # Error handling
  if (length(mu)!=d){
    stop("Dimension mismatch: mu and Sigma!")
  }
  
  Schol <- chol(Sigma$S)
  tmp <- matrix(rnorm(n*p),n,p)
  tmp <- tmp%*%Schol
  tmp <- tmp%*%t(Sigma$L) + apply(t(Sigma$N),2,FUN=function(x){
    sqrt(x)*rnorm(n)
  })
  t(apply(tmp,1,FUN=function(x){
    x+mu
  }))
}
