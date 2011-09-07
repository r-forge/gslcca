\name{readSpectra}
\alias{readSpectra}
\title{
Read in Power Spectra and Auxiliary Information
}
\description{
Read in a tab-delimited text file containing power spectra and auxiliary information, i.e. subject, treatment and time variables.
}
\usage{
readSpectra(file, info = 1:3, treatment = NULL, ...)
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
  \item{\dots}{
    Further arguments passed to \code{\link{read.delim}}
}
}
\details{
Information variables that are read as character vectors are converted to 
factors, using the levels passed to the \code{treatment} argument for the 
treatment variable (if found). Note these levels must match the original values 
of the treatment variable, this argument simply specifies the ordering of the levels.
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
\keyword{ IO }
\keyword{ files }
