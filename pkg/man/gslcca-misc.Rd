\name{gslcca-misc}
\alias{signatures}
\alias{fitted.gslcca}
\title{
Functions to Access Components of the Results of a GSLCCA Analysis
}
\description{
These functions provide ways to access key components of the \code{"gslcca"}
object returned by the function \code{gslcca}.
}
\usage{
signatures(object, ...)

\method{fitted}{gslcca}(object, ...)
}
\arguments{
  \item{object}{
    a \code{"gslcca"} object as returned by the function \code{gslcca}.
}
  \item{\dots}{
  arguments to pass on
}
}
\value{
\code{signatures} returns the signatures estimated in a GSLCCA analysis. These 
are the coefficients or loadings of the multivariate response. If a subject 
factor was specified in the original analysis, the function returns a data frame 
with one column per subject, otherwise a single vector is returned. 

The \code{fitted} method returns the (projected) fitted values from the GSLCCA,
i.e. the x scores.
}
\author{
Heather Turner
}
\seealso{
 \code{\link{gslcca}}
}
\examples{
data(clonidine)
banded <- bandSpectra(clonidine$spectra, breaks = c(4, 8, 13, 30), 
             labels = c("delta", "theta", "alpha", "beta", "gamma"))
result <- gslcca(banded, "Critical Exponential",
    time = Time, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine, subset = Rat == "42")
    
signatures(result)    
head(fitted(result))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
\keyword{ misc }% __ONLY ONE__ keyword per line
