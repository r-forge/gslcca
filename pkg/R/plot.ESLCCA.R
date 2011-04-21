plot.ESLCCA <- function(x, type = "signature", individual = TRUE, overlay = FALSE,
                        ask = dev.interactive(), lattice = FALSE,
                        main = NULL, xlab = NULL, ylab = NULL,
                        col = NULL, lty = NULL, pch = 1, legend.x = "topright",
                        space = "bottom", corner = NULL, columns = 2, ...){
    ## control over legend and title?

    ns <- nlevels(x$subject)
    ls <- levels(x$subject)
    nt <- nlevels(x$treatment)
    lt <- levels(x$treatment)

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
    if (fitted)
        dat <- data.frame(subject = x$subject, trt = x$treatment, time = x$time, y = x$y, fit = x$fitted.values)

    if (missing(xlab)) xlab <- ifelse(signature, "Frequency", "Time")
    if (missing(ylab)) ylab <- ifelse(signature, "Loading", "Response")
    if (missing(col)) col <- rainbow_hcl(ifelse(signature, ns, nt))
    if (missing(lty)) lty <- seq_len(ifelse(signature, ns, nt))


    if (individual) {
        if (overlay) {
            if (signature) {
                matplot(x$frequencies,t(x$signatures),col=col,type='l',
                        xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Signatures corresponding to different subjects', main))
                legend(x = legend.x,y=NULL, legend = paste('Subject',ls), col = col,
                       lty = lty, pch = NULL, merge = TRUE)
            }
            else if (fitted) { #show fitted curves only else plot too noisy
                fit <- cast(dat, time ~ subject + trt, value = "fit")
                matplot(fit[, 1], fit[,-1], col = col, type = "l", lty,
                        xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Fitted values corresponding to different subjects', main))
                legend(x = legend.x,y=NULL, legend = lt, col = col,
                       lty = lty, pch = NULL, merge = TRUE)
            }
        }
        else if (lattice) {
            ## need to explicitly print plots as not typing directly into console
            if (signature) {
                sig <- t(x$signatures)
                print(xyplot(sig ~ c(row(sig)) | colnames(sig)[col(sig)], type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Signatures corresponding to different subjects', main),
                             xlab= xlab, ylab = ylab, as.table = TRUE, ...))
            }
            else if (fitted) {
                ## plotting of both points and lines not compatible with group
                print(xyplot(x$fitted.values ~ x$time | x$subject, group=x$treatment, type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Fitted values corresponding to different subjects', main),
                             xlab = xlab,  ylab = ylab, as.table = TRUE,
                             key = list(space = space, corner = corner, lines = list(col = col, lty = lty),
                             text = list(levels(x$treatment)), border = TRUE, columns = columns)), ...)
                nr <- round(sqrt(ns))
                for (i in seq_len(ns)) {
                    id <- x$subject == ls[i]
                    trellis.focus("panel", column = (i - 1) %% nr + 1, row = (i - 1) %/% nr + 1)
                    panel.points(x$time[id], x$y[id], col = col, cex = 0.8, pch = pch)
                    trellis.unfocus()
                }
            }
        }
        else {
            for (i in seq_len(ns)) {
                if (signature) {
                    plot(x$frequencies,x$signatures[i,],xlab=xlab,ylab=ylab,type='l')
                    title(ifelse(missing(main), paste('Signature for subject',ls[i]),
                                 paste(main,ls[i])))
                }
                else if (fitted) {
                    obs <- cast(dat, time ~ trt, subset =  x$subject == ls[i], value = "y")
                    fit <- cast(dat, time ~ trt, subset =  x$subject == ls[i], value = "fit")
                    matplot(obs[, 1], obs[,-1], col = col, type = "p", pch = pch,
                            xlab = xlab,  ylab = xlab, ...)
                    matlines(fit[, 1], fit[, -1], col = col, type = "l", lty = lty)
                    title(ifelse(missing(main), paste('Fitted values for subject',ls[i]),
                                 paste(main,ls[i])))
                    legend(x = legend.x, y=NULL, legend = lt, col = col,
                           lty = lty, pch = NULL, merge = TRUE)
                }
            }
        }
    }
    else {
        if (signature) {
            plot(colMeans(x$signatures),xlab=xlab, ylab = ylab, type='l', ...)
            title(ifelse(missing(main), 'Mean Signature', main))
        }
        else if (fitted) {
            fit <- cast(dat, time ~ trt, value = "fit", fun.aggregate = mean)
            matplot(fit[, 1], fit[,-1], col = col, type = "l", lty = lty,
                    xlab = xlab,  ylab = ylab, ...)
            title(ifelse(missing(main), 'Mean fitted curves', main))
            legend(x = legend.x,y=NULL, legend = lt, col = col, lty = lty,
                   pch = NULL, merge = TRUE)
       }
    }
}
