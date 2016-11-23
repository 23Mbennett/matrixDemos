\name{RREF}
\Rdversion{1.1}
\alias{RREF}
\title{
Reduced Row-Echelon Form.
}

\description{
This function uses Gaussian elimination to 
transform a matrix to reduced row-echelon form. The actual computations are
done by the \code{\link{GaussianElimination}} function.
}

\usage{
RREF(X, ...)
}

\arguments{
  \item{X}{a numeric matrix.}
  \item{\dots}{arguments to be passed to \code{\link{GaussianElimination}}.}
}


\value{
the reduced row-echelon form of \code{X}.
}

\references{
John Fox (2009)
\emph{A Mathematical Primer for Social Statistics.}
Thousand Oaks CA: Sage.
}

\author{
John Fox
}

\seealso{
\code{\link{GaussianElimination}}
}

\examples{
A <- matrix(c(1,2,3,4,5,6,7,8,10), 3, 3) # a nonsingular matrix
A
RREF(A) # the reduced row-echelon form of A

b <- 1:3
RREF(A, b)  # solving the matrix equation Ax = b

RREF(A, diag(3)) # inverting A

B <- matrix(1:9, 3, 3) # a singular matrix
B
RREF(B)

RREF(B, b)

RREF(B, diag(3))
}
\keyword{algebra}
