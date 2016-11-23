cholesky <- function(X, tol=sqrt(.Machine$double.eps)){
    # returns the Cholesky square root of the nonsingular, symmetric matrix X
    # tol: tolerance for checking for 0 pivot
    # algorithm from Kennedy & Gentle (1980)
    if (!is.numeric(X)) stop("argument is not numeric")
    if (!is.matrix(X)) stop("argument is not a matrix")
    n <- nrow(X)
    if (ncol(X) != n) stop("matrix is not square")
    if (max(abs(X - t(X))) > tol) stop("matrix is not symmetric")
    D <- rep(0, n)
    L <- diag(n)
    i <- 2:n
    D[1] <- X[1, 1]
    if (abs(D[1]) < tol) stop("matrix is numerically singular")
    L[i, 1] <- X[i, 1]/D[1]
    for (j in 2:(n - 1)){
        k <- 1:(j - 1)
        D[j] <- X[j, j] - sum((L[j, k]^2) * D[k])
        if (abs(D[j]) < tol) stop("matrix is numerically singular")
        i <- (j + 1):n
        L[i, j] <- (X[i, j] -
                        colSums(L[j, k] * t(L[i, k, drop=FALSE]) * D[k]))/D[j]
        }
    k <- 1:(n - 1)
    D[n] <- X[n, n] - sum((L[n, k]^2) * D[k])
    if (abs(D[n]) < tol) stop("matrix is numerically singular")
    L %*% diag(sqrt(D))
    }

