\name{bandSpectra}
\alias{bandSpectra}
\title{
Average Spectra within Frequency Bands
}
\description{
Given a matrix of spectra, where each row gives the power at each frequency, 
create banded spectra by averaging across columns within frequency bands.
}
\usage{
bandSpectra(spectra, breaks = NULL, labels = NULL, ...)
}
\arguments{
  \item{spectra}{
    A matrix of spectra.
}
  \item{breaks}{
    The break points between the frequency bands, including a lower limit
    (exclusive) and an upper limit (inclusive).
}
  \item{labels}{
    Labels for the frequency bands.
}
  \item{\dots}{
    Further arguments passed to cut.
}
}
\value{
A matrix of banded spectra.
}
\author{
Heather Turner
}
\seealso{
\code{\link{gslcca}}
}
\examples{
## apply to arbitrary matrix to average over column groups
M <- matrix(1:12, 3, 4)
colnames(M) <- 1:4
bandSpectra(M, breaks = c(0, 2, 4))

## apply to clonidine experiment
data(clonidine)

banded <- bandSpectra(clonidine$spectra, breaks = c(4, 8, 10, 13, 15, 20, 28), 
             labels = c("delta", "theta", "alpha1", "alpha2", "not named",
                        "beta1", "beta2", "gamma"))
## run gslcca
result <- gslcca(banded, "Critical Exponential",
    time = Time, treatment = Treatment, subject = Rat, separate = FALSE,            
    subject.smooth = 4, data = clonidine)
                 
plot(result, "fitted", lattice = TRUE)
plot(result, "signatures", overlay = TRUE)
}
\keyword{ manip }
