\name{clonidine}
\alias{clonidine}
\docType{data}
\title{
  EEG Experiment Investigating the Effect of Clonidine in Rats
}
\description{
Power spectra derived from electroencephalogram (EEG) measurements of
brain activity for eight rats from a balanced four-way cross-over design
experiment in which rats were dosed with three doses of clonidine (low,
medium, or high) and vehicle.
}
\usage{data(clonidine)}
\format{
  A data frame consisting of the following
  \describe{
    \item{Rat}{ a factor with eight levels: "35", "36", ..., "42". }
    \item{Treatment}{ a factor with four levels: "Control", "Low Dose",
      "Middle Dose", "High Dose". }
    \item{Time}{ a variable giving the time in seconds at the start of
      the five minute time interval over which the mean spectra was
      calculated. }
    \item{spectra}{ a matrix whose rows contain mean power spectra corresponding
    to the design points specified by the first three variables. }
  }
}
\details{
Eight rats were implanted with radiotelemetry transmitters
intraperitoneally for recording via electroencephalogram (EEG). After a
recovery period, the rats were orally dosed with three different doses of
clonidine and vehicle, each rat receiving all four dose levels in a
balanced four way 4 crossover design.

The rats were singly housed on a standard 12:12 h light-dark cycle and
received standard diet and water ad libitum. The rats were dosed at the
begining of the light period and these data correspond to following 12
hours. At least 48 hours were allowed before the next treatment was
administered to avoid any drug carryover effect. 

The difference in the voltages between the two electrodes at each time
point was carried forward for further analysis. This consisted of taking
12-second time slices and subjecting them to FFT, which gave a power
spectrum for each time slice over 36 frequencies. The power spectra were
filtered by removing spectra where power at any frequency was greater
than 0.1 powers greater than this were regarded as artefacts). The
remaining spectra were averaged over five minute time intervals (so up
to 25 spectra were averaged) to give the final set of mean spectra.

To reduce the size of the data set distributed with the package, the
data have been subsampled, so that every other mean spectra has been
kept. 
}
\source{
Original data produced by Pfizer Global R & D, UK.
}
\seealso{
    \code{\link{readSpectra}} for creating data structures of this type from 
    tab-delimited data files.
    
    \code{\link{gslcca}} for GSLCCA analysis.
}
\examples{
data(clonidine)

######################
## Example periodogram
######################

eg <- with(clonidine, { 
           list(Rat == "42" & Time == 0 & Treatment == "Control",
                Rat == "42" & Time == 8400 & Treatment == "Control")
           })

freq <- as.numeric(colnames(clonidine$spectra))
plot(freq, clonidine$spectra[eg[[2]]], type = "l", ylab =
    expression(paste("Power (", mu, V^2, ")")), lty = 2,
    xlab = "Frequency (Hz)", main = "Examples of frequency periodograms")
lines(1:36, clonidine$spectra[eg[[1]]], lty = 1)
legend("topright", c("0-300 seconds", "8400-8700 seconds"), lty = 1:2)

#######################################################################
## GSLCCA analysis for Rat 42 using separate Critical Exponential model
## for each treatment
#######################################################################

### Unsmoothed data
result <- gslcca(spectra, "Critical Exponential",
    time = Time, treatment = Treatment, subject.smooth = FALSE,
    data = clonidine, subset = Rat == "42")

## drug signature
plot(result, "signature")

## projected values
plot(result, "fitted")

### Smoothed data - automatically select number of roots
result <- gslcca(spectra, "Critical Exponential",
    time = Time, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine, subset = Rat == "42")

## drug signature
plot(result, "signature")

## projected values
plot(result, "fitted")

## projected values
plot(result, "scores")

}
\keyword{datasets}
