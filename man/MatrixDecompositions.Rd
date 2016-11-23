\name{MatrixDecompositions}
\alias{MatrixDecompositions}
\alias{eig}
\alias{QR}
\alias{SVD}
\alias{cholesky}

\title{
Various matrix decompositions.
}
\description{
\code{cholesky} computes the Cholesky square root of a square, symmetric, nonsingular matrix;
\code{eig} computes the eigenvalues and eigenvectors of a square, symmetric matrix; 
\code{QR} computers the QR decomposition of a matrix; 
and \code{SVD} computes the singular-value decomposition of a matrix
}
\usage{
cholesky(X, tol = sqrt(.Machine$double.eps))

eig(X, tol = sqrt(.Machine$double.eps), max.iter = 100, retain.zeroes = TRUE)

QR(X, tol = sqrt(.Machine$double.eps))

SVD(X, tol = sqrt(.Machine$double.eps))
}

\arguments{
  \item{X}{for \code{cholesky}, a square, symmetric, nonsingular matrix; 
    for \code{eig}, a square symmetric matrix; for \code{QR} or \code{SVD}, a matrix.}
  \item{tol}{0 tolerance; defaults to the square-root of the machine double-epsilon}
  \item{max.iter}{iteration limit.}
  \item{retain.zeroes}{retain 0 eigenvalues and corresponding eigenvectors? Defaults to \code{TRUE}.}
}
\details{
The Cholesky square root of a square symmetric nonsingular matrix is found using the method of Kennedy and Gentle (1980). 
The singular-value decomposition of a matrix $X$ is computed from the eigenstructure of $X^{T}X$. 
The eigenvalues and vectors computed by \code{eig} are found via the iterated QR decomposition. 
The QR decomposition is computed by Graham-Schmidt orthonormalization. 
}

\value{
\code{cholesky} returns the Cholesky square root matrix; 
\code{eig} returns a list with elements \code{values} and \code{vectors}; 
\code{QR} returns a list with elements \code{Q}, \code{R}, and \code{rank}; 
\code{SVD} returns a list with elements \code{d} (the singular values), \code{U} and \code{V}.
}

\author{
John Fox and Georges Monette
}

\references{
J. Fox (2009)
\emph{A Mathematical Primer for Social Statistics.}
Thousand Oaks CA: Sage.

W. J. Kennedy Jr. and J. E. Gentle (1980)
\emph{Statistical Computing.}
New York: Dekker.
}

%% note: the LaTeX comment charater % must be escaped in the examples

\examples{
(C <- matrix(c(1,2,3,2,5,6,3,6,10), 3, 3)) # nonsingular, symmetric
cholesky(C)
cholesky(C) \%*\% t(cholesky(C))  # check

(EC <- eig(C)) # eigenanalysis of C
EC$vectors \%*\% diag(EC$values) \%*\% t(EC$vectors) # check

(D <- matrix(1:6, 3, 2))

(SD <- SVD(D))  # singular value decomposition of D
SD$U \%*\% diag(SD$d) \%*\% t(SD$V) # check

(QRD <- QR(D))  # QR decomposition of D
QRD$Q \%*\% QRD$R  # check

# principal-components analysis
(R <- cor(DavisThin))
eig(R)
}

\keyword{algebra}

