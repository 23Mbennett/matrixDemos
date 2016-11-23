\name{GaussianElimination}
\Rdversion{1.1}
\alias{GaussianElimination}

\title{
Gaussian Elimination
}
\description{
This function transforms its first matrix argument to reduced row-echelon form by
Gaussian elimination, applying the same row operations to its second matrix argument
(if it is present).
The function can be used for a variety of purposes, such as inverting square matrices, and
solving square or rectangular systems of linear simultaneous equations.
}
\usage{
GaussianElimination(A, B, tol = sqrt(.Machine$double.eps), 
    verbose = FALSE, fractions = FALSE)
}

\arguments{
  \item{A}{numeric matrix to be transformed to reduced row-echelon form.}
  \item{B}{optional numeric matrix to be brought along for the ride.}
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
a matrix giving the reduced row-echelon form of \code{A} adjoined to the corresponding
transformed \code{B}.
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
\code{\link[MASS]{fractions}}
}

\examples{
A <- matrix(c(1,2,3,4,5,6,7,8,10), 3, 3) # a nonsingular matrix
A
GaussianElimination(A, diag(3)) # find inverse of A

b <- 1:3
GaussianElimination(A, b, verbose=TRUE, fractions=TRUE)  # solve the matrix equation Ax = b
}
\keyword{algebra}
