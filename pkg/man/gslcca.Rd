\name{gslcca}
\alias{gslcca}
\title{ Perform a Generalised Semi Linear Canonical Correlation Analysis }
\description{
  Performs (Generalised) Semi Linear Canonical Correlation Analysis,
  i.e. computes the canonical correlation between a data matrix and a 
  nonlinear function of time. SLCCA is extended by allowing parameters
  to vary by a treatment factor and allowing adjustment for covariates.
  Shortcuts are provided for PK/PD models suitable for analysing EEG data.
}
\usage{
gslcca(Y, formula = "Double Exponential", time, subject = NULL,
    treatment = NULL, ref = 1, separate = TRUE, partial = ~1, data, subset,
    global.smooth = FALSE, subject.smooth = TRUE, pct.explained = 0.98,
    start = NULL, ...)
}
\arguments{
  \item{Y}{ a data matrix. }
  \item{formula}{ either a nonlinear model formula specifying a function
    of \code{time} including parameters, or one of the character strings
    \code{"Double Exponential"} or \code{"Critical Exponential"} to specify the
    corresponding PK/PD model (see Details). }
  \item{time}{ a vector of time values corresponding to the rows of \code{Y}. }
  \item{subject}{ an optional factor grouping the rows of Y. If
    specified, data corresponding to to each level will be analysed
    separately (after any smoothing specified by
    \code{global.smooth}). }
  \item{treatment}{ an optional factor (nested within \code{subject} if
    specified) identifying groups for which separate parameters
    should be estimated. }
  \item{ref}{ the reference level of \code{treatment} for
    which \code{formula} is set to \code{~ 0} (may be \code{NULL}). }
  \item{separate}{ if \code{treatment} is specified, whether to estimate
    separate parameters of \code{formula} and separate canonical coefficients
    for each treatment or just separate canonical coefficients. }
  \item{partial}{ a linear model formula specifying covariates to
    partial out of the CCA analysis (may be \code{NULL}).}
  \item{data}{ an optional data frame in which to evaluate the variables
    in \code{time}, \code{subject}, \code{treatment} and \code{partial}. }
  \item{subset}{}
  \item{global.smooth}{ controls the smoothing of \code{Y} via SVD. Must
    be one of \code{FALSE} for no smoothing, \code{TRUE} for smoothing
    with an automatically selected number of roots, or a scalar specifying the
    number of roots to use. }
  \item{subject.smooth}{ controls the smoothing of \code{Y} within the
    levels of \code{subject}. Accepts the same values as
    \code{global.smooth}. }
  \item{pct.explained}{ a scalar between 0 and 1 indicating the desired
    minimum percentage of variance explained by the SVD approximation
    when \code{subject.smooth} or \code{global.smooth} is \code{TRUE}.}    
  \item{start}{ a named list of starting values for the parameters. Each
    element must have the same length as the number of non-reference
    levels of \code{treatment} or have length one, in which case the
    same starting values are used for each non-reference treatment
    level. May be \code{NULL} when \code{formula} is a character string,
    however it is recommended that starting values are always
    provided.}
  \item{...}{ arguments passed on to \code{optim}.}
}
\details{
  The function fits the following model:
  \deqn{\mathbf{Y(t)} \times \mathbf{A} = \mathbf{X(t, \theta)} \times
    \mathbf{B}}{
    Y(t) * A = X(t, theta) * B}
  where \eqn{\mathbf{Y(t)}}{Y(t)} is a data
  matrix with rows of observations recorded at times
  \eqn{\mathbf{t}}{t}, \eqn{\mathbf{A}}{A} is a vector of loadings,
  \eqn{\mathbf{X(t, \theta)}}{X(t, theta)} is a matrix with columns
  containing a nonlinear function with unknown parameters for each
  non-reference treatment level, and \eqn{\mathbf{B}}{B} is a
  vector of coefficients.

  The parameters \eqn{\mathbf{A}, \theta}{A, theta} and
  \eqn{\mathbf{B}}{B} are estimated to optimise the correlation
  between the left- and right-hand sides of the model.

  If \code{partial} specifies a matrix of covariates, \eqn{G}, to be
  partialled out, then the canonical correlation analysis is based on
  the residuals from the multivariate linear models
  \code{lm(Y ~ 0 + G)}
  and
  \code{lm(X ~ 0 + G)}
  When \code{partial = ~1} this is equivalent to centering the columns
  of \eqn{Y} and \eqn{X}.

  The nonlinear function defining \eqn{X} may be specified by a
  character string for two pharmacokinetic/pharmacodynamic models:
  \describe{
    \item{\code{"Double Exponential"}}{\code{~ exp(-time/exp(K1))-exp(-time/exp(K2))}}
    \item{\code{"Critical Exponential"}}{\code{~ time*exp(-time/exp(K1))}}
  }

  The data matrix \code{Y} may be pre-smoothed to remove global
  artefacts. This is achieved by approximating the data from the
  number of SVD roots specified by \code{global.smooth}. If the
  number of roots is not specified explicitly, \eqn{k} roots are
  selected such that 
  \eqn{\lambda_k/\lambda_1 >0.001}{lambda_k/lambda_1 >0.001} and
  \eqn{\frac{\sum_{j=1}^k \lambda_j}{\sum_{j=1}^r \lambda_j} \ge p}{
    (sum_{j=1}^k lambda_j)/(sum_{j=1}^r lambda_j) >= p}
   where \eqn{k=1,\ldots,r}{k = 1, ..., r}, \eqn{r} is the rank of
  \code{Y}, \eqn{\lambda_1, \lambda_2, \ldots}{lambda_1, lambda_2,
    ...} are the eigenvalues of the variance matrix and \eqn{p} is
  given by \code{pct.explained}.

  The data are then analysed for each level of \code{subject} separately
  and the data may again be smoothed at this level to remove local
  artefacts as specified by \code{local.smooth}.

  The coefficients of the left-hand side matrix are sometimes referred
  to as \dfn{loadings} and are described in Brain et al (2011) as
  \dfn{signatures}.

  The projected values \eqn{\mathbf{Y(t)} \times \mathbf{A}}{Y(t) * A}
  and \eqn{\mathbf{X(t, \theta)} \times \mathbf{B}}{X(t, theta) * B} are
  known as the y scores and x scores respectively. In this context, they
  can be thought of as the (projected) observed and fitted values.
  
  Signatures, observed and fitted values can be visualised using
  \code{\link{plot.gslcca}}. 
}
\value{
  An object of class \code{"gslcca"}, which is a list with components
  \item{call}{ the call to \code{gslcca}.  }
  \item{ycoef}{ a data frame giving the coefficients of the left-hand
    side matrix for each subject.  }
  \item{xcoef}{  a data frame giving the coefficients of the right-hand
    side matrix for each subject. }
  \item{scores}{ a data frame containing the variables specified by
  \code{subject}, \code{treatment} and \code{time} along with new variables
  \code{yscores} and \code{xscores}, containing the projected values
  for all subjects.} 
  \item{ref}{ the reference level of \code{treatment}.}
  \item{nonlinear.parameters}{ a data frame giving the estimated
    nonlinear parameters for each subject. }
  \item{y}{ a list giving the left-hand side matrix for each subject. }
  \item{x}{ a list giving the right-hand side matrix for each subject. }
  \item{global.smooth}{ the value of \code{global.smooth}. }
  \item{subject.smooth}{ the value of \code{subject.smooth}. }
  \item{opt}{ the output from \code{optim} for each subject. }
}
\references{
  Brain, P., Strimenopoulou, F. and Ivarsson, M. (2011). Analysing
  electroencephalogram (EEG) data using Extended Semi-Linear Canonical
  Correlation Analysis.  Submitted. }
\author{ Foteini Strimenopoulou and Heather Turner }
\seealso{ \code{\link{plot.gslcca}}, \code{\link{varySmooth}}}
\examples{
data(clonidine)
}
\keyword{ multivariate }
