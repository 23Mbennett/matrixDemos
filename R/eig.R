eig <- function(X, tol=sqrt(.Machine$double.eps), max.iter=100, retain.zeroes=TRUE){
  # returns the eigenvalues and eigenvectors of a square, symmetric matrix using the iterated QR decomposition
  # X: a square, symmetric matrix
  # tol: 0 tolerance
  # max.iter: iteration limit
  # retain.zeroes: retain 0 eigenvalues?
  if (!is.numeric(X) || !is.matrix(X) || nrow(X) != ncol(X) || any(abs(X - t(X)) > tol)) 
      stop("X must be a numeric, square, symmetric matrix")
  i <- 1
  Q <- diag(nrow(X))
  while (i <= max.iter){
    qr <- QR(X, tol=tol)
    Q <- Q %*% qr$Q 
    X <- qr$R %*% qr$Q
    if (max(abs(X[lower.tri(X)])) <= tol) break
    i <- i + 1
  }
  if (i > max.iter) warning("eigenvalues did not converge")
  values <- diag(X)
  if (!retain.zeroes){
    nonzero <- values != 0
    values <- values[nonzero]
    Q <- Q[, nonzero, drop=FALSE]
  }
  list(values=values, vectors=Q)
}
