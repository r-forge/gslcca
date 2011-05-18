\name{plot.ESLCCA}
\alias{plot.ESLCCA}
\title{
  Plot Results from an ESLCCA Analysis}
}
\description{
  Plot \dQuote{signatures} (loadings of the Y matrix) or (projected)
  observed and fitted values (Y and X scores) from an \dQuote{ESLCCA} object.
}
\usage{
    plot.ESLCCA(x, type = "signature", individual = TRUE,
        overlay = FALSE, ask = dev.interactive(), lattice = FALSE,
        main = NULL, xlab = NULL, ylab = NULL, col = NULL, lty = NULL,
        pch = NULL, legend.x = "topright", space = "bottom",
        corner = NULL, columns = 2, ...)
}
\arguments{
  \item{x}{ an object of class \code{"ESLCCA"}. }
  \item{type}{ the type of plot: either \code{"signature"} or \code{"fitted"}. }
  \item{individual}{ logical: if \code{TRUE} values are plotted for
    each level of \code{x$subject} individually, otherwise values are
    averaged over subjects.  }
  \item{overlay}{ logical:  if \code{TRUE} the values for each subject
    are overlaid on the same plot, otherwise a separate plot is
    produced for each subject. Ignored when \code{individual} is \code{FALSE}.}
  \item{ask}{ logical: if \code{TRUE}, the user is \emph{ask}ed before each plot, see
    \code{\link{par}(ask=.)}. Ignored when only one plot is produced. }
  \item{lattice}{ logical: if \code{TRUE} the values for each subject
    are displayed in separate panels of a trellis plot. }
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{legend.x}{
%%     ~~Describe \code{legend.x} here~~
}
  \item{space}{
%%     ~~Describe \code{space} here~~
}
  \item{corner}{
%%     ~~Describe \code{corner} here~~
}
  \item{columns}{
%%     ~~Describe \code{columns} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, type = "signature", individual = TRUE, overlay = FALSE,
                        ask = dev.interactive(), lattice = FALSE,
                        main = NULL, xlab = NULL, ylab = NULL,
                        col = NULL, lty = NULL, pch = NULL, legend.x = "topright",
                        space = "bottom", corner = NULL, columns = 2, ...){
    ## control over legend and title?

    ns <- nlevels(x$scores$subject)
    ls <- levels(x$scores$subject)
    nt <- nlevels(x$scores$treatment)
    lt <- levels(x$scores$treatment)

    ## count plots and set ask to FALSE if not needed (also reset on exit if necessary)
    nplots <- length(type)
    if (nplots > 1) stop("Only one type of plot may be plotted at a time")
    if (individual && !(overlay | lattice)) nplots <- ns
    if (nplots <= prod(par("mfcol"))) ask <- FALSE
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    signature <- length(agrep("signature", type))
    fitted <- length(agrep("fitted", type))
    if (!(signature | fitted))
        stop("\"type\" not recognised - must be \"signature\" or \"fitted\"")
    if (signature)
        freq <- seq_len(ncol(x$ycoef))

    if (missing(xlab)) xlab <- ifelse(signature, "Frequency", "Time")
    if (missing(ylab)) ylab <- ifelse(signature, "Coefficient", "Score")
    if (missing(col)) col <- c(hex(polarLUV(H = 0, C = 0, L = 60)), rainbow_hcl(ifelse(signature, ns - 1, nt - 1), c = 100, l = 60)) # assumes ref
    if (missing(lty)) lty <- seq_len(ifelse(signature, ns, nt))
    if (missing(pch)) pch <- seq_len(ifelse(signature, ns, nt))


    if (individual) {
        if (overlay) {
            if (signature) {
                matplot(freq, t(x$ycoef), col=col,type='l', xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Signatures corresponding to different subjects', main))
                legend(x = legend.x,y=NULL, legend = paste('Subject',ls), col = col,
                       lty = lty, pch = NULL, merge = TRUE)
            }
            else if (fitted) { #show fitted curves only else plot too noisy
                fit <- cast(x$scores, time ~ subject + treatment, value = "xscores")
                matplot(fit[, 1], fit[,-1], col = col, type = "l", lty = lty, pch = pch,
                        xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Fitted values corresponding to different subjects', main))
                legend(x = legend.x,y=NULL, legend = lt, col = col,
                       lty = lty, pch = NULL, merge = TRUE)
            }
        }
        else if (lattice) {
            ## need to explicitly print plots as not typing directly into console
            if (signature) {
                sig <- t(x$ycoef)
                print(xyplot(sig ~ c(row(sig)) | rownames(x$ycoef)[col(sig)], type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Signatures corresponding to different subjects', main),
                             xlab= xlab, ylab = ylab, as.table = TRUE, ...))
            }
            else if (fitted) {
                ## plotting of both points and lines not compatible with group
                print(xyplot(xscores ~ time | subject, group=treatment, type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Fitted values corresponding to different subjects', main),
                             xlab = xlab,  ylab = ylab, as.table = TRUE, data = x$scores,
                             key = list(space = space, corner = corner, lines = list(col = col, lty = lty,  pch = pch),
                             type = "b", text = list(lt), border = TRUE, columns = columns), ...))
                layout <- trellis.currentLayout()
                nc <- ncol(layout)
                for (i in seq_len(ns)) {
                    trellis.focus("panel", column = (i - 1) \%\% nc + 1, row = (i - 1) \%/\% nc + 1)
                    with(x$scores, {
                        id <- subject == ls[i]
                        panel.points(time[id], yscores[id], col = col[treatment[id]], pch = pch[treatment[id]], cex = 0.6)
                        })
                    trellis.unfocus()
                }
            }
        }
        else {
            for (i in seq_len(ns)) {
                if (signature) {
                    plot(freq, x$ycoef[i,],xlab=xlab,ylab=ylab,type='l')
                    title(ifelse(missing(main), paste('Signature for subject',ls[i]),
                                 paste(main,ls[i])))
                }
                else if (fitted) {
                    yscores <- cast(x$scores, time ~ treatment, subset =  subject == ls[i], value = "yscores")
                    xscores <- cast(x$scores, time ~ treatment, subset =  subject == ls[i], value = "xscores")
                    matplot(yscores[, 1], yscores[,-1], col = col, type = "p", pch = pch,
                            xlab = xlab,  ylab = xlab, cex = 0.6, ...)
                    matlines(xscores[, 1], xscores[, -1], col = col, type = "l", lty = lty)
                    title(ifelse(missing(main), paste('Fitted values for subject',ls[i]),
                                 paste(main,ls[i])))
                    legend(x = legend.x, y=NULL, legend = lt, col = col,
                           lty = lty, pch = pch, merge = TRUE)
                }
            }
        }
    }
    else {
        if (signature) {
            plot(colMeans(x$ycoef),xlab=xlab, ylab = ylab, type='l', ...)
            title(ifelse(missing(main), 'Mean Signature', main))
        }
        else if (fitted) {
            xscores <- cast(x$scores, time ~ treatment, value = "xscores", fun.aggregate = mean)
            matplot(xscores[, 1], xscores[,-1], col = col, type = "l", lty = lty,
                    xlab = xlab,  ylab = ylab, ...)
            title(ifelse(missing(main), 'Mean fitted curves', main))
            legend(x = legend.x,y=NULL, legend = lt, col = col, lty = lty,
                   pch = NULL, merge = TRUE)
       }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
