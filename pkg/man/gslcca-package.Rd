\name{gslcca-package}
\alias{gslcca-package}
\docType{package}
\title{ Generalised Semi-Linear Canonical Correlation Analysis }
\description{
GSLCCA estimates the parameters of a given nonlinear model to maximize
the correlation with a linear combination of multiple response
variables.  
}
\details{
  GSLCCA was developed to characterize brain activity monitored using EEG
  (electroencephalogram) under different treatment regimens in research
  and development projects at Pfizer.
  
  EEG results are typically summarised as a matrix of power spectra,
  giving the power (microvolts squared) over multiple frequencies (Hz)
  recorded at different times. GSLCCA is then used to find a linear
  combination of the powers that is maximally correlated to a nonlinear
  pharmacometric model.
  
  The \code{gslcca} package can be used to relate any multivariate
  response to any function of time, however it provides shortcuts for
  specifying standard PK/PD models (Critical Exponential or Double
  Exponential). It also has options to analyse subjects separately;
  to estimate parameters separately for different treatments and to
  partial out covariates prior to analysis.
  
  Utility functions are provided to plot the results and to investigating
  the effect of pre-smoothing on the analysis.
}
\author{
Foteini Strimenopoulou and Heather Turner

Maintainer: Heather Turner <ht@heatherturner.net>
}
\references{
  Brain, P., Strimenopoulou, F. and Ivarsson, M. (2011). Analysing
  electroencephalogram (EEG) data using Extended Semi-Linear Canonical
  Correlation Analysis.  Submitted.
}

\keyword{ package }
\seealso{
\code{\link{gslcca}} for the main function, with links to associated functions.
}
