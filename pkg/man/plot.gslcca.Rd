\name{plot.gslcca}
\alias{plot.gslcca}
\title{
  Plot Results from a GSLCCA Analysis
}
\description{
  Plot \sQuote{signatures} (loadings of the \eqn{\mathbf{Y}}{Y} matrix) or (projected)
  observed and fitted values (\eqn{\mathbf{Y}}{Y} and
  \eqn{\mathbf{X}}{X} scores) from an \dQuote{gslcca} object. 
}
\usage{
    \method{plot}{gslcca}(x, type = "signature", series = x$treatment,
        mean = FALSE, overlay = (type == "signature"), ask = dev.interactive(),
        lattice = (type != "signature"), main = NULL, xlab = NULL, ylab = NULL,
        col = NULL, lty = 1, lwd = 1.5, pch = NULL, legend.x = "topright",
        space = "bottom", corner = NULL, columns = 2, ...)
}
\arguments{
  \item{x}{ an object of class \code{"gslcca"}. }
  \item{type}{ the type of plot: either \code{"signature"},
    \code{"fitted"} or \code{"scores"}. }
  \item{series}{ if \code{x$treatment} is \code{NULL} and there are
    multiple time series for each subject, a factor to differentiate
    these series. }
  \item{mean}{ logical: if \code{TRUE} values are 
    averaged over subjects.  }
  \item{overlay}{ for \code{type = \"signature\"} only, a logical value
    specifying whether to overlay signatures for different subjects on
    the same plot or create a separate plot for each subject. }
  \item{ask}{ logical: if \code{TRUE}, the user is \emph{ask}ed before each plot, see
    \code{\link{par}(ask=.)}. Ignored when only one plot is produced. }
  \item{lattice}{ logical: if \code{TRUE} the values for each subject
    are displayed in separate panels of a trellis plot. }
  \item{main}{ an overall title for the plot. }
  \item{xlab}{ a title for the x axis. }
  \item{ylab}{ a title for the y axis. }
  \item{col}{ a vector specifying a colour palette to use for the
    lines/points. By default the palette is specified by
    \code{\link{rainbow_hcl}} if available, otherwise the standard
    plotting colours are used. }
  \item{lty}{ a vector of line types, see \code{\link{par}}.}
  \item{lwd}{ a vector of line widths, see \code{\link{par}}.}
  \item{pch}{ a vector of plotting characters, see \code{\link{points}}.}
  \item{legend.x}{ a character vector to use in the legend.}
  \item{space}{ when \code{type = "fitted"} and \code{lattice = TRUE},
    the space outside the plot region in which to put the legend, see
    \code{\link{xyplot}}. Overrides the \code{corner} argument. } 
  \item{corner}{ when \code{type = "fitted"} and \code{lattice = TRUE},
    the corner inside the plot region in which to put the  legend, see
    \code{\link{xyplot}}. }  
  \item{columns}{ when \code{type = "fitted"} and \code{lattice = TRUE},
    how many columns to use to format the legend, see \code{\link{xyplot}}. }
  \item{\dots}{ arguments passed on to the plotting functions
    (\code{\link{plot}} for single curves, \code{\link{xyplot}} for
    trellis graphics, else \code{\link{matplot}})}
}
\details{
  For \code{type = "signature"}, the loadings of the \eqn{\mathbf{Y}}{Y}
  matrix from the GSLCCA analysis are displayed for each subject or
  their average is displayed, as specified by \code{individual}. When both
  \code{individual} and \code{overlay} are \code{TRUE}, the arguments
  \code{col}, \code{lty} and \code{legend.x} control the line style and
  legend text for each subject.

  For \code{type = "fitted"}, the \eqn{\mathbf{X}}{X} scores are
  displayed with lines. If there are multiple \eqn{\mathbf{X}}{X} scores
  for each time point, then a line is displayed for each times series
  specified by \code{series} (usually the treatment factor used in
  GSLCCA) and the arguments \code{col}, \code{lty},
  \code{pch} and \code{legend.x} control the display for each
  series. By default, grey or black is used to display values
  corresponding to the reference treatment (specified by \code{x$ref}). 
  
  For \code{type = "scores"} the \eqn{\mathbf{X}}{X} scores are displayed as 
  for \code{type = "fitted"} and points are added for the \eqn{\mathbf{Y}}{Y}
  scores, again coloured by \code{series}.
}
\author{
  Heather Turner
}
\seealso{
\code{\link{gslcca}}
}
\examples{
data(clonidine)

### Separate critical exponential model for each treatment
result <- gslcca(spectra, "Critical Exponential",
    time = Time, treatment = Treatment, subject = Rat,
    subject.smooth = 4, data = clonidine)

## plot of individual signatures
plot(result, type = "signature", overlay = FALSE)

## plot of individual projected observed and fitted values
plot(result, type = "scores", lattice = FALSE)

## plot of individual projected fitted values
plot(result, type = "fitted", lattice = FALSE)

## plot of individual signatures - overlaid
plot(result, type = "signature")

## plot of individual signatures - lattice
plot(result, type = "signature", overlay = FALSE, lattice = TRUE)

## plot of individual projected observed and fitted values - lattice
plot(result, type = "scores")

## plot of individual projected fitted values - lattice
plot(result, type = "fitted")

## plot of mean signature
plot(result, type = "signature", mean = TRUE)

## plot of mean fitted
plot(result, type = "fitted", mean = TRUE)

}
\keyword{ hplot }