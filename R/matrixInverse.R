matrixInverse <- function(X, tol=sqrt(.Machine$double.eps), ...){
    # returns the inverse of nonsingular X
    if ((!is.matrix(X)) || (nrow(X) != ncol(X)) || (!is.numeric(X))) 
        stop("X must be a square numeric matrix")
    n <- nrow(X)
    X <- GaussianElimination(X, diag(n), tol=tol, ...) # append identity matrix
        # check for 0 rows in the RREF of X:
    if (any(apply(abs(X[,1:n]) <= sqrt(.Machine$double.eps), 1, all)))
        stop ("X is numerically singular")
    X[,(n + 1):(2*n)]  # return inverse
    }

