QR <- function(X, tol=sqrt(.Machine$double.eps)){
  # QR decomposition by Graham-Schmidt orthonormalization
  # X: a matrix
  # tol: 0 tolerance
  if (!is.numeric(X) || !is.matrix(X)) stop("X must be a numeric matrix")
  length <- function(u) sqrt(sum(u^2))
  U <- X
  E <- matrix(0, nrow(X), ncol(X))
  E[, 1] <- U[, 1]/length(U[, 1])
  for (j in 2:ncol(U)){
    for (k in 1:(j - 1)){
      U[, j] <- U[, j] - (X[, j] %*% E[, k]) * E[, k]
    }
    len.U.j <- length(U[, j])
    if (len.U.j > tol) E[, j] <- U[, j]/len.U.j
  }
  R <- t(E) %*% X
  R[abs(R) < tol] <- 0
  rank <- sum(rowSums(abs(R)) > 0)
  list(Q=-E, R=-R, rank=rank) # negated to match qr()
}
