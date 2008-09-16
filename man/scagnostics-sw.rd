\name{scagnostics}
\alias{scagnostics}
\alias{scagnostics.default}
\alias{scagnostics.matrix}
\alias{scagnostics.data.frame}
\alias{scagnostics_2d}
\title{Calculcate scagnostics for a scatterplot}
\author{Hadley Wickham <h.wickham@gmail.com>}

\description{
Scagnostics summarise potentially interesting patterns in 2d scatterplot
}
\usage{scagnostics(x, ...)}
\arguments{
\item{x}{object to calculate scagnostics on: a vector, a matrix or a data.frame}
\item{...}{...}
}

\details{Current scagnostics are:

\itemize{
\item Outlying
\item Skewed
\item Clumpy
\item Sparse
\item Striated
\item Convex
\item Skinny
\item Stringy
\item Monotonic
}

These are described in more detail in:
Graph-Theoretic Scagnostics, Leland Wilkinson, Anushka
Anand, Robert Grossman.
\url{http://www.ncdm.uic.edu/publications/files/proc-094.pdf}

You can call the function with two 1d vectors to get a single vector
of scagnostics, or with a 2d structure (matrix or data frame) to get
scagnostics for every combination of the variables.}

\examples{scagnostics(1:10, 1:10)
scagnostics(rnorm(100), rnorm(100))
scagnostics(mtcars)
scagnostics(as.matrix(mtcars))

if (require(rggobi)) ggobi(scagnostics(mtcars))}
\keyword{hplot}
