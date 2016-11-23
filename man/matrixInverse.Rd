\name{matrixInverse}
\Rdversion{1.1}
\alias{matrixInverse}

\title{Matrix Inverse}

\description{
This function computes the inverse of a square nonsingular matrix
by the method of Gaussian elimination. The actual computations are
performed by the \code{\link{GaussianElimination}} function.
}

\usage{
matrixInverse(X, tol = sqrt(.Machine$double.eps), ...)
}

\arguments{
  \item{X}{a square numeric matrix}
  \item{tol}{A small tolerance to be used to detect numbers that are
  effectively 0; the default is the square-root of machine precision for
  double-precision floating-point numbers.}
  \item{\dots}{aguments to be passed down to \code{\link{GaussianElimination}}.}
}

\value{
the matrix inverse of \code{X}, if it exists.
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
matrixInverse(A, fractions=TRUE)
round(matrixInverse(A) \%*\% A, 6) # check
}

\keyword{algebra}
