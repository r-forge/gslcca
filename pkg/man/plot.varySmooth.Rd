\name{plot.varySmooth}
\alias{plot.varySmooth}
\title{
Plot Results from GSLCCA Run with Different Smoothing Settings
}
\description{
Plot the optimised value, \sQuote{signatures} (loadings of the
\eqn{\mathbf{Y}}{Y} matrix) or (projected) observed and fitted values
(\eqn{\mathbf{Y}}{Y} and \eqn{\mathbf{X}}{X} scores) for a
\code{varySmooth} object.
}
\usage{
\method{plot}{varySmooth}(x, type = "opt", series = x[[1]]$treatment,
  subject = levels(x[[1]]$subject), ask = dev.interactive(),
  main = NULL, xlab = NULL, ylab = NULL, col = NULL, lty = 1, lwd = 1.5, 
  pch = NULL, space = "bottom", corner = NULL, columns = 2, ...)
}
\arguments{
  \item{x}{ an object of class \code{"varySmooth"} as returned by \code{\link{varySmooth}}}
  \item{type}{ the type of plot: either \code{"opt"}, \code{"signature"} or \code{"fitted"}. }
  \item{series}{ if \code{x[[1]]$treatment} is \code{NULL} and there are
    multiple time series for each subject, a factor to differentiate
    these series. }
  \item{subject}{ the subjects to include in the plot, specified by
    levels of \code{x$subject}. }
  \item{ask}{logical: if \code{TRUE}, the user is \emph{ask}ed before each plot, see
    \code{\link{par}(ask=.)}. Ignored when only one plot is produced.}
  \item{main}{ an overall title for the plot.}
  \item{xlab}{ a title for the x axis.}
  \item{ylab}{ a title for the y axis. }
  \item{col}{ a vector specifying a colour palette to use for the
    lines/points. By default the palette is specified by
    \code{\link{rainbow_hcl}} if available, otherwise the standard
    plotting colours are used.}
  \item{lty}{ a vector of line types, see \code{\link{par}}.}
  \item{lwd}{ a vector of line widths, see \code{\link{par}}.}
  \item{pch}{a vector of plotting characters, see \code{\link{points}}.}
  \item{space}{ when \code{type = "fitted"},
    the space outside the plot region in which to put the legend, see
    \code{\link{xyplot}}. Overrides the \code{corner} argument.}
  \item{corner}{ when \code{type = "fitted"},
    the corner inside the plot region in which to put the  legend, see
    \code{\link{xyplot}}.}
  \item{columns}{when \code{type = "fitted"},
    how many columns to use to format the legend, see \code{\link{xyplot}}. }
  \item{\dots}{ further arguments passed to \code{\link{xyplot}}. }
}
\details{
  For \code{type = "opt"}, the value of the optimisation criterion,
  \eqn{\log(1 - R^2)}{log(1 - R^2)} is plotted against the number of
  roots used in the GSLCCA, trellised by subject.

  For \code{type = "signature"}, the loadings of the \eqn{\mathbf{Y}}{Y}
  matrix from the GSLCCA are displayed for each subject,
  trellised by the number of roots used in the GSLCCA. The arguments
  \code{col} and \code{lty} control the line style and
  legend text for each subject.

  For \code{type = "fitted"}, the \eqn{\mathbf{X}}{X} scores are
  displayed with lines and the \eqn{\mathbf{Y}}{Y} scores are displayed
  with points and the display is again trellised by the number of roots
  used in the GSLCCA. If there are multiple \eqn{\mathbf{X}}{X} scores
  for each time point, then a line is displayed for each times series
  specified by \code{series} (usually the treatment factor used in
  GSLCCA) and the arguments \code{col}, \code{lty} and \code{pch}
  control the display for each series. By default, grey or black is used
  to display values corresponding to the reference treatment (specified
  by \code{x$ref}). 
}
\author{ Heather Turner }
\seealso{
\code{\link{gslcca}}, \code{\link{varySmooth}}
}
\examples{}
\keyword{ multivariate }