Ginv <- function(A, tol=sqrt(.Machine$double.eps), verbose=FALSE, 
        fractions=FALSE){
    # return an arbitrary generalized inverse of the matrix A
    # A: a matrix
    # tol: tolerance for checking for 0 pivot
    # verbose: if TRUE, print intermediate steps
    # fractions: try to express nonintegers as rational numbers
    if (fractions) {
        mass <- requireNamespace("MASS", quietly=TRUE)
        if (!mass) stop("fractions=TRUE needs MASS package")
    }
    m <- nrow(A)
    n <- ncol(A)
    B <- GaussianElimination(A, diag(m), tol=tol, verbose=verbose, 
        fractions=fractions)
    L <- B[,-(1:n)]
    AR <- B[,1:n]
    C <- GaussianElimination(t(AR), diag(n), tol=tol, verbose=verbose, 
        fractions=fractions)
    R <- t(C[,-(1:m)])
    AC <- t(C[,1:m])
    ginv <- R %*% t(AC) %*% L
    if (fractions) MASS::fractions (ginv) else round(ginv, round(abs(log(tol, 10))))
    }

