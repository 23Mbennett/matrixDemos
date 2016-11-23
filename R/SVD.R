SVD <- function(X, tol=sqrt(.Machine$double.eps)){
  # compute the singular-value decomposition of a matrix X from the eigenstructure of X'X
  # X: a matrix
  # tol: 0 tolerance
  VV <- eig(t(X) %*% X, tol=tol, retain.zeroes=FALSE)
  V <- VV$vectors
  d <- sqrt(VV$values)
  U <- X %*% V %*% diag(1/d,nrow=length(d)) # magically orthogonal
  list(d=d, U=U, V=V)
}
