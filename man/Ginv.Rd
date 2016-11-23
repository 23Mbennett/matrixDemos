\name{Ginv}
\Rdversion{1.1}
\alias{Ginv}

\title{
Generalized Matrix Inverse.
}

\description{
Computes a generalized matrix inverse of a square or rectangular matrix by
Gaussian elimination. The actual computations are
performed by the \code{\link{GaussianElimination}} function.
}

\usage{
Ginv(A, tol = sqrt(.Machine$double.eps), verbose = FALSE, fractions = FALSE)
}
\arguments{
  \item{A}{a square or rectangular numeric matrix.}
  \item{tol}{small tolerance to detect numbers that are to be considered 0;
  defaults to the machine predicion for double-precision floating-point
  numbers.}
  \item{verbose}{if \code{TRUE} (the default is \code{FALSE}), show the
  progress of the Gaussian elimination algorith.}
  \item{fractions}{if \code{TRUE} (the default is \code{FALSE}), try
  to convert non-integers to rational fractions, using the 
  \code{\link[MASS]{fractions}} function in the \pkg{MASS} package.}
}


\value{
a matrix giving a generalized inverse of \code{A}.
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

B <- matrix(1:9, 3, 3) # a singular matrix

Ginv(A, fractions=TRUE)  # a generalized inverse of A = inverse of A
round(Ginv(A) \%*\% A, 6)  # check

Ginv(B, fractions=TRUE)  # a generalized inverse of B
B \%*\% Ginv(B) \%*\% B   # check

C <- matrix(1:12, 3, 4) # rectangular matrix
Ginv(C, fractions=TRUE)
round(C \%*\% Ginv(C) \%*\% C, 5)   # check
}

\keyword{algebra}
