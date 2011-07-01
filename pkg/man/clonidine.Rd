\name{clonidine}
\alias{clonidine}
\docType{data}
\title{
  EEG Experiment Investigating the Effect of Clonidine in Rats
}
\description{
Power spectra derived from electroencephalogram (EEG) measurements of
brain activity for two rats from a balanced four-way cross-over design
experiment in which rats were dosed with three doses of clonidine (low,
medium, or high) and vehicle.
}
\usage{data(clonidine)}
\format{
  A list of two components
  \item{spectra} a matrix of mean power spectra, see Details.
  \item{design} a data frame of variables specifying the design point
  corresponding to each spectrum
  \describe{
    \item{Rat} a factor with two levels: "41" and "42".
    \item{Treatment} a factor with four levels: "Control", "Low Dose",
    "Middle Dose", "High Dose".
    \item{From} a variable giving the time in seconds at the start of
    the five minute time interval over which the mean spectra was
    calculated. 
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
}
\source{
Original data produced by Pfizer Global R \& D, UK.
}
\examples{
data(clonidine)

######################
## Example periodogram
######################

with(clonidine$design, {
    eg1 <-  Rat == "42" & From == 0 & Treatment == "Control"
    eg2 <-  Rat == "42" & From == 8400 & Treatment == "Control"
})

## Subject 64 control dose, Time 1 = 0, Time 2 = 8400
plot(1:36, clonidine$spectra[eg2], type = "l", ylab = "Power", lty = 2,
    xlab = "Frequency", main = "Examples of frequency periodograms")
lines(1:36, clonidine$spectra[eg1], lty = 1)
legend("topright", c("0-300 seconds", "8400-8700 seconds"), lty = 1:2)

#######################################################################
## ESLCCA analysis for Rat 42 using separate Critical Exponential model
## for each treatment
#######################################################################

### Unsmoothed data
result <- ESLCCA(clonidine$spectra, "Critical Exponential",
    time = From, treatment = Treatment, subject.smooth = FALSE,
    data = clonidine$design, subset = Rat == "42")

## drug signature
plot(result, "signature")

## projected values
plot(result, "fitted")

### Smoothed data - automatically select number of roots
result <- ESLCCA(clonidine$spectra, "Critical Exponential",
    time = From, treatment = Treatment, subject.smooth = TRUE,
    data = clonidine$design, subset = Rat == "42")

## drug signature
plot(result, "signature")

## projected values
plot(result, "fitted")

### Vary number of roots and plot optimised value
result <- list()
for (r in 2:15) {
    result[[r]] <- ESLCCA(clonidine$spectra, "Critical Exponential",
        time = From, treatment = Treatment, subject.smooth = r,
        data = clonidine$design, subset = Rat == "42")
}

result <- result[-1]

## possibly 7 roots may be needed here
opt <- c(sapply(result, function(x) sapply(x$opt, "[[", "value")))
xyplot(opt ~ 2:15, type="b",
       main = "Log(1 - Correlation)",
       xlab= "Roots", ylab = "Log(1 - Corr)", as.table = TRUE)
}

xscores <- unlist(lapply(result, "[[", "xscores"))
yscores <- unlist(lapply(result, "[[", "yscores"))
time <- unlist(lapply(result, "[[", "time"))
treatment <- unlist(lapply(result, "[[", "treatment"))
roots <- rep(2:15, each = length(result[[1]]$xscores))

xyplot(xscores ~ time | roots, group = treatment, type = "l", col = 1:4, lty = 1:4,
    xlab = "Roots",  ylab = "Score", as.table = TRUE)
layout <- trellis.currentLayout()
nc <- ncol(layout)
for (i in 1:14) {
  trellis.focus("panel", column = (i - 1) %% nc + 1, row = (i - 1) %/% nc + 1)
  id <- roots == i + 1
  panel.points(time[id], yscores[id], col = (1:4)[treatment[id]],
  pch = (1:4)[treatment[id]], cex = 0.6)
  trellis.unfocus()
}


\keyword{datasets}
