GaussianElimination <- function(A, B, tol=sqrt(.Machine$double.eps),
                                verbose=FALSE, fractions=FALSE){
  # A: coefficient matrix
  # B: right-hand side vector or matrix
  # tol: tolerance for checking for 0 pivot
  # verbose: if TRUE, print intermediate steps
  # fractions: try to express nonintegers as rational numbers
  # If B is absent returns the reduced row-echelon form of A.
  # If B is present, reduces A to RREF carrying B along.
  printMatrix <- function(A){
    if (fractions) print(fractions(as.matrix(A)))
    else print(round(as.matrix(A), round(abs(log(tol, 10)))))
  }
  frac <- function(x) as.character(fractions(x))
  rowadd <- function(x, from, to, mult) {
    x[to, ] <- x[to, ] + mult * x[from, ]
    x
  }
  rowswap <- function(x, from, to) {
    x[c(to, from), ] <- x[c(from, to), ]
    x
  }
  rowmult <- function(x, row, mult) {
    x[row, ] <- mult * x[row, ]
    x
  }
  if ((!is.matrix(A)) || (!is.numeric(A)))
    stop("argument must be a numeric matrix")
  n <- nrow(A)
  m <- ncol(A)
  if (!missing(B)){
    B <- as.matrix(B)
    if (!(nrow(B) == nrow(A)) || !is.numeric(B))
      stop("argument must be numeric and must match the number of row of A")
    A <- cbind(A, B)
  }
  i <- j <- 1
  if (verbose){
    cat("\nInitial matrix:\n")
    printMatrix(A)
  }
  while (i <= n && j <= m){
    if (verbose) cat("\nrow:", i, "\n")
    while (j <= m){
      currentColumn <- A[,j]
      currentColumn[1:n < i] <- 0
      # find maximum pivot in current column at or below current row
      which <- which.max(abs(currentColumn))
      pivot <- currentColumn[which]
      if (abs(pivot) <= tol) { # check for 0 pivot
        j <- j + 1
        next
      }
      if (which > i) {
        A <- rowswap(A, i, which) # exchange rows (E3)
        if (verbose) {
          cat("\n exchange rows", i, "and", which, "\n")
          printMatrix(A)
        }
      }
      A <- rowmult(A, i, 1/pivot) # pivot (E1)
      if (verbose && abs(pivot - 1) > tol){
        cat("\n multiply row", i, "by", 
            if (fractions) frac(1/pivot) else 1/pivot, "\n")
        printMatrix(A)
      }
      for (k in 1:n){
        if (k == i) next
        factor <- A[k, j]
        if (abs(factor) < tol) next
        A <- rowadd(A, i, k, -factor) # sweep column j (E2)
        if (verbose){
          if (abs(factor - 1) > tol){
            cat("\n multiply row", i, "by",
                if (fractions) frac(abs(factor)) else abs(factor),
                if (factor > 0) "and subtract from row" else "and add to row", k, "\n")
          }
          else{
            if (factor > 0) cat("\n subtract row", i, "from row", k, "\n")
            else cat("\n add row", i, "from row", k, "\n")
          }
          printMatrix(A)
        }
      }
      j <- j + 1
      break
    }
    i <- i + 1
  }
  # 0 rows to bottom
  zeros <- which(apply(A[,1:m], 1, function(x) max(abs(x)) <= tol))
  if (length(zeros) > 0){
    zeroRows <- A[zeros,]
    A <- A[-zeros,]
    A <- rbind(A, zeroRows)
  }
  rownames(A) <- NULL
  ret <- if (fractions) fractions(A) else round(A, round(abs(log(tol, 10))))
  if (verbose) invisible(ret) else ret
}
