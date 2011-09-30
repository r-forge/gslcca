\name{snapshot}
\alias{snapshot}
\alias{plot.snapshot}
\title{
Perform Snapshot Analysis Comparing Two Sets of Signatures from GSLCCA
}
\description{
Performs t-tests to test for differences in the mean loadings between two sets 
of signatures found using GSLCCA. Corrects for multiple comparisons.
}
\usage{
snapshot(x, y, p.adjust.method = "fdr", normalise = TRUE, ...)
}
\arguments{
  \item{x, y}{
    a matrix of signatures from a GSLCCA analysis, with one column per subject.
}
  \item{normalise}{
    logical value indicating whether or not to normalise the signatures first,
    so that the sum of squares is equal to one.
}
  \item{p.adjust.method}{
    the method of adjustment for multiple comparisons, see 
    \code{\link{p.adjust}}.
}
  \item{\dots}{
    further arguments passed to \code{rowttests}
}
}
\details{
In GSLCCA, a signature is the set of loadings or coefficients, \eqn{\mathbf{a}}, 
of the multivariate response \eqn{\mathbf{Y}}, that give the highest correlation 
between the y scores \eqn{\mathbf{Ya}} and the fitted nonlinear model.

In the context of EEG analysis the signature represents the relative importance
of each frequency in relation to the PK/PD model.

In a GSLCCA analysis with multiple subjects, a signature will have been estimated
for each subject. This function compares the mean loading of each variable 
between two such sets of signatures, using t-tests.

When \code{normalise=TRUE}, each signature is first scaled so that the sum of
squared loadings is equal to one.

The p-values from the t-test are corrected for multiplicity, using the method 
given by \code{p.adjust.method}, which specifies the false 
discovery rate by default (Benjamini & Hochberg, 1995).

}
\value{
A list with the following components
\item{call}{ the call to \code{snapshot}.}
\item{x }{ the first set of (normalised) signatures. }
\item{y }{ the second set of (normalised) signatures. }
\item{pvalue }{ the corrected pvalue. }
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a 
practical and powerful approach to multiple testing. 
\emph{Journal of the Royal Statistical Society Series B}, \bold{57}, 289--300.
}
\author{
Heather Turner
}
\seealso{
\code{\link{gslcca}}, \code{\link{plot.gslcca}}
}
\examples{
data(clonidine)

### Arbitrarily split the rats into two groups of four
### - expect results to be roughly the same!
grp1 <- subset(clonidine, Rat < 39)
grp2 <- subset(clonidine, Rat >= 39)

## fit same model to each group
result <- gslcca(spectra, "Double Exponential",
    time = Time, subject = Rat, treatment = Treatment, 
    subject.smooth = TRUE, pct.explained = 0.96, data = grp1)

result2 <- gslcca(spectra, "Double Exponential",
    time = Time, subject = Rat, treatment = Treatment, 
    subject.smooth = TRUE, pct.explained = 0.96, data = grp2)

## run snapshot analysis and plot 
## - observed differences far from significant as expected
snap <- snapshot(signatures(result), signatures(result2))
par(mfrow = c(2,2))
plot(snap, type = c("signatures", "means", "pvalue"), 
     names = c("Group 1", "Group 2"))

### Artificially create a difference by swapping the data
### from the delta and theta bands for Group 1
grp1$spectra <- grp1$spectra[, c(5:8, 1:4, 9:36)]
colnames(grp1$spectra) <- colnames(grp2$spectra)

## refit model for group 1 and re-run snapshot analysis
result <- gslcca(spectra, "Double Exponential",
    time = Time, subject = Rat, treatment = Treatment, 
    subject.smooth = TRUE, pct.explained = 0.96, data = grp1)    
snap <- snapshot(signatures(result), signatures(result2))

## now shows highly significant difference in delta,
## significant difference in theta
## (Group 1 signatures less consistent here)
par(mfrow = c(2,2))
plot(snap, type = c("signatures", "means", "pvalue"), 
     names = c("Group 1", "Group 2"))

## compact display
par(mfrow = c(1,1))
plot(snap, type = "compact", names = c("Group 1", "Group 2"))    

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
