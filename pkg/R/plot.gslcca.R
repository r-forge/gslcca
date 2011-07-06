plot.gslcca <- function(x, type = "signature", series = x$treatment, individual = TRUE,
                        overlay = FALSE, ask = dev.interactive(), lattice = FALSE,
                        main = NULL, xlab = NULL, ylab = NULL,
                        col = NULL, lty = NULL, pch = NULL, legend.x = "topright",
                        space = "bottom", corner = NULL, columns = 2, ...){
    ## control over legend and title?
    subject <- x$subject
    one <- FALSE
    if (is.null(subject)) {
        lattice <- overlay <- FALSE
        individual <- one <- TRUE
        subject <- factor(rep.int(1, length(x$xscores)))
    }
    ns <- nlevels(subject)
    ls <- levels(subject)

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
        freq <- seq_len(nrow(x$ycoef))

    treatment <- series
    nt <- nlevels(treatment)
    lt <- levels(treatment)
    if (is.null(treatment)) {
        ## create factor to id series
        time <- as.numeric(paste(unclass(subject), x$time, sep = ""))
        ord <- order(time)
        reps <- table(time)
        treatment <- factor(sequence(reps)[order(ord)])
        nt <- 1
    }

    if (missing(xlab)) xlab <- ifelse(signature, "Frequency", "Time")
    if (missing(ylab)) ylab <- ifelse(signature, "Coefficient", "Score")

    if (missing(col)) {
        if (exists("rainbow_hcl", mode= "function"))
            col <- c(hex(polarLUV(H = 0, C = 0, L = 60)), rainbow_hcl(ifelse(signature, ns - 1, nt - 1), c = 100, l = 60)) # assumes ref
        else
            col <- seq_len(ifelse(signature, ns, nt))
    }
    if (missing(lty)) lty <- seq_len(ifelse(signature, ns, nt))
    if (missing(pch)) pch <- seq_len(ifelse(signature, ns, nt))


    if (individual) {
        if (overlay) {
            if (signature) {
                matplot(freq, x$ycoef, col=col,type='l', xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Signatures Corresponding to Different Subjects', main))
                legend(x = legend.x,y=NULL, legend = paste('Subject',ls), col = col,
                       lty = lty, pch = NULL, merge = TRUE)
            }
            else if (fitted) { #show fitted curves only else plot too noisy
                fit <- tapply(x$xscores, list(x$time, interaction(subject, treatment)), I)
                matplot(names(fit[, 1]), fit[,-1], col = col, type = "l", lty = lty, pch = pch,
                        xlab= xlab, ylab= ylab, ...)
                title(ifelse(missing(main), 'Fitted Values Corresponding to Different Subjects', main))
                if (nt > 1) {
                    legend(x = legend.x,y=NULL, legend = lt, col = col,
                           lty = lty, pch = NULL, merge = TRUE)
                }
            }
        }
        else if (lattice) {
            ## need to explicitly print plots as not typing directly into console
            if (signature) {
                sig <- as.matrix(x$ycoef)
                print(xyplot(sig ~ c(row(sig)) | colnames(sig)[col(sig)], type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Signatures Corresponding to Different Subjects', main),
                             xlab= xlab, ylab = ylab, as.table = TRUE, ...))
            }
            else if (fitted) {
                ## plotting of both points and lines not compatible with group
                if (nt > 1) {
                    key <- list(space = space, corner = corner, lines = list(col = col, lty = lty,  pch = pch),
                                type = "b", text = list(lt), border = TRUE, columns = columns)
                }
                else {
                    key <- NULL
                    pch <- rep(pch[1], nlevels(treatment))
                    col <- rep(col[1], nlevels(treatment))
                }
                print(xyplot(x$xscores ~ x$time | subject, group=treatment, type=c("l"), col = col, lty = lty,
                             main = ifelse(missing(main), 'Fitted Values Corresponding to Different Subjects', main),
                             xlab = xlab,  ylab = ylab, as.table = TRUE, key = key, ...))
                layout <- trellis.currentLayout()
                nc <- ncol(layout)
                for (i in seq_len(ns)) {
                    trellis.focus("panel", column = (i - 1) %% nc + 1, row = (i - 1) %/% nc + 1)
                    id1 <- subject == ls[i]
					id2 <- unclass(as.factor(treatment))[id1]
                    panel.points(x$time[id1], x$yscores[id1], col = col[id2],
                                 pch = pch[id2], cex = 0.6)
                    trellis.unfocus()
                }
            }
        }
        else {
            if (signature) {
                for (i in seq_len(ns)) {
                    plot(freq, x$ycoef[,i],xlab=xlab,ylab=ylab,type='l', col = col[1], lty = lty[1], ...)
                    title(ifelse(missing(main), ifelse(one, 'Signature', paste('Signature for Subject',ls[i])),
                                 paste(main,ls[i])))
                }
            }
            if (fitted) {
                yscores <- tapply(x$yscores, list(x$time, treatment, subject), I)
                xscores <- tapply(x$xscores, list(x$time, treatment, subject), I)
                for (i in seq_len(ns)) {
                    matplot(names(yscores[, 1, i]), yscores[, , i], col = col, type = "p", pch = pch,
                            xlab = xlab,  ylab = xlab, cex = 0.6, ...)
                    matlines(names(xscores[, 1, i]), xscores[, , i], col = col, type = "l", lty = lty)
                    title(ifelse(missing(main), ifelse(one, 'Fitted Values', paste('Fitted Values for Subject',ls[i])),
                                 paste(main,ls[i])))
                    if (nt > 1) {
                        legend(x = legend.x, y=NULL, legend = lt, col = col,
                               lty = lty, pch = pch, merge = TRUE)
                    }
                }
            }
        }
    }
    else {
        if (signature) {
            plot(rowMeans(x$ycoef),xlab=xlab, ylab = ylab, type='l', col = col[1], lty = lty[1], ...)
            title(ifelse(missing(main), 'Mean Signature', main))
        }
        else if (fitted) {
            xscores <- tapply(x$xscores, list(x$time, treatment), I)
            matplot(rownames(xscores), xscores, col = col, type = "l", lty = lty,
                    xlab = xlab,  ylab = ylab, ...)
            title(ifelse(missing(main), 'Mean Fitted Values', main))
            if (nt > 1) {
                legend(x = legend.x,y=NULL, legend = lt, col = col, lty = lty,
                       pch = NULL, merge = TRUE)
            }
       }
    }
    invisible()
}
