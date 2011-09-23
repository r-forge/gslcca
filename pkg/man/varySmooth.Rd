\name{varySmooth}
\alias{varySmooth}
\title{
Vary the Level of Smoothing used Prior to GSLCCA
}
\description{
Re-run GSLCCA with different amounts of smoothing at the subject level
prior to analysis.
}
\usage{
varySmooth(x, subject.smooth = 1:10, ...)
}
\arguments{
  \item{x}{ a \code{"gslcca"} object as returned by \code{\link{gslcca}}.}
  \item{subject.smooth}{ a vector of integers specifying the number of
    SVD roots use when smoothing the data}
  \item{\dots}{ currently ignored}
}
\details{
This function can be used to investigate the effect of smoothing on
GSLCCA.

When the \code{subject.smooth} argument of \code{\link{gslcca}} is an
integer, the data matrix for each subject is approximated using the
corresponding number of SVD roots.
}
\value{
  An object of class \code{"varySmooth"} which is a list of
  \code{"gslcca"} objects obtained by re-running the original GSLCCA
  given by \code{x} with each value of \code{subject.smooth}
}
\author{ Heather Turner }
\seealso{
\code{\link{gslcca}}, \code{\link{plot.varySmooth}}
}
\examples{
data(clonidine)

### Smoothed data - automatically select number of roots
result <- gslcca(spectra, "Critical Exponential",
    time = Time, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine, subset = Rat == "42")

### Vary number of roots 
multiRoots <- varySmooth(result, 2:15)

## plot optimised value
plot(multiRoots, "opt")

## plot fitted values
plot(multiRoots, "fitted")

## plot signature
plot(multiRoots, "signature")
}
\keyword{ multivariate }
