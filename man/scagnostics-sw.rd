\name{scagnostics}
\alias{scagnostics}
\alias{scagnostics.default}
\alias{scagnostics.data.frame}
\title{Calculcate scagnostics for a scatterplot}
\author{Hadley Wickham <h.wickham@gmail.com>}

\description{
Scagnostics summarise potentially interesting patterns in 2d scatterplot
}
\usage{scagnostics(x,...)}
\arguments{
\item{x}{object to calculate scagnostics on}
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
\url{http://www.ncdm.uic.edu/publications/files/proc-094.pdf}}

\examples{scagnostics(1:10, 1:10)
scagnostics(rnorm(100), rnorm(100))
scagnostics(mtcars)

if (require(rggobi)) ggobi(scagnostics(mtcars))}
\keyword{hplot}
