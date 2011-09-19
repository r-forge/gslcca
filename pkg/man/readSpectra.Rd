\name{readSpectra}
\alias{readSpectra}
\title{
Read in Power Spectra and Auxiliary Information
}
\description{
Read in a tab-delimited text file containing power spectra and auxiliary information, i.e. subject, treatment and time variables.
}
\usage{
readSpectra(file, info = 1:3, 
    treatment = c("Control", "Low Dose", "Middle Dose", "High Dose"), 
    resolution = 1, end = 42900, nfreq = 36, ...)
}
\arguments{
  \item{file}{
    The name of the file.
}
  \item{info}{
    Numeric indices of columns containing auxiliary information - the remaining
    columns are assumed to contain the power spectra.
}
  \item{treatment}{
    A vector of levels for the treatment variable, assumed to be named 
    "treatment" (case insensitive) - ignored if no such variable is found.
}
  \item{resolution}{
    The desired resolution of frequencies, in Hz. Data read in at a higher 
    resolution is aggregated.  
}
  \item{end}{
    The last time point to be included in the result - observations at further 
    time points are dropped from the data.
  }
  \item{nfreq}{
    The number of frequencies to be included in the result (after any 
    aggregation specified by \code{resolution}), i.e. the number of columns in 
    the returned matrix of spectra.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{read.delim}}
}
}
\details{
Information variables that are read as character vectors are converted to 
factors, using the levels passed to the \code{treatment} argument for the 
treatment variable (if found). Note these levels must match the original values 
of the treatment variable, this argument simply specifies the ordering of the
levels.
}
\value{
A data frame containing the information variables with their original names and 
the matrix of spectra, named "spectra".
}
\author{
Heather Turner
}
\seealso{
\code{\link{gslcca}}, \code{\link{bandSpectra}}
}
\examples{
\dontrun{
## Read data in from tab-delimited text file
DATA <- readSpectra("Compound1.txt", info = 1:3, 
    treatment = c("Control", "Low Dose", "Middle Dose", "High Dose"),
    resolution = 1, end = 42900, nfreq = 36)

## Use result directly in call to gslcca
result <- gslcca(spectra,  "Double Exponential",
                 time=Time, subject = Rat, 
                 treatment = Treatment, data = DATA,
                 separate = FALSE, global = FALSE,
                 subject.smooth = TRUE, pct.explained = 0.96)
}
}
\keyword{ IO }
\keyword{ files }
