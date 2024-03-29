plot.gslcca <- function(x, type = "signature", series = x$treatment,
                        mean = FALSE, 
                        overlay = length(agrep(type, "signature")),
                        ask = dev.interactive(), 
                        lattice = !length(agrep(type, "signature")), 
                        main = NULL, xlab = NULL, ylab = NULL, 
                        col = NULL, lty = 1, lwd = 1.5,
                        pch = NULL, legend.x = "topright", space = "bottom",
                        corner = NULL, columns = 2, ...){
    
    ## control over legend and title?
    subject <- x$subject
    one <- FALSE
    if (is.null(subject)) {
        mean <- lattice <- overlay <- FALSE
        one <- TRUE
        subject <- factor(rep.int(1, length(x$xscores)))
    }
    ns <- nlevels(subject)
    ls <- levels(subject)

    ## count plots and set ask to FALSE if not needed (reset on exit if necessary)
    nplots <- length(type)
    if (nplots > 1) stop("Only one type of plot may be plotted at a time")
    if (!mean && !(overlay | lattice)) nplots <- ns
    if (nplots <= prod(par("mfcol"))) ask <- FALSE
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    signature <- length(agrep("signature", type))
    scores <- length(agrep("scores", type))
    fitted <- length(agrep("fitted", type))
    if (!(signature | scores | fitted))
        stop("\"type\" not recognised - must be \"signature\", \"scores\" or \"fitted\"")
    if (signature) {
        xaxis <- list()
        xaxis$freq.nm <- rownames(x$ycoef)
        freq <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", xaxis$freq.nm)))
        if(any(is.na(freq))) freq <- xaxis$at <- seq_along(freq)
        else xaxis <- NULL
    }

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

    if (missing(xlab)) xlab <- ifelse(signature, "Frequency (Hz)", "Time")
    if (missing(ylab)) ylab <- ifelse(signature, "Coefficient", "Score")

    if (missing(col)) {
        if (exists("rainbow_hcl", mode= "function"))
            col <- c(hex(polarLUV(H = 0, C = 0, L = 60)),
                     rainbow_hcl(ifelse(signature, ns - 1, nt - 1),
                                 c = 100, l = 60)) # assumes ref
        else
            col <- seq_len(ifelse(signature, ns, nt))
    }
    if (missing(pch)) pch <- seq_len(ifelse(signature, ns, nt))


    if (!mean) {
        if (overlay) {
            if (signature) {
                matplot(freq, x$ycoef, col=col,type='l', xlab= xlab, ylab= ylab,
                        lty = lty, lwd = lwd,  xaxt = "n", ...)
                axis(1, xaxis$at, xaxis$freq.nm)
                title(ifelse(missing(main),
                             'Signatures Corresponding to Different Subjects',
                             main))
                legend(x = legend.x,y=NULL, legend = paste('Subject',ls), col = col,
                       lty = lty, lwd = lwd, pch = NULL, merge = TRUE)
            }
            else stop("overlay only implemented for type = \"signature\"")
        }
        else if (lattice) {
            ## need to explicitly print plots as not typing directly into console
            if (signature) {
                sig <- as.matrix(x$ycoef)
                print(xyplot(sig ~ c(row(sig)) | colnames(sig)[col(sig)],
                             type=c("l"), col = col, lty = lty, lwd = lwd,
                             scales = list(x = xaxis),
                             main = ifelse(missing(main),
                                           'Signatures Corresponding to Different Subjects', main),
                             xlab= xlab, ylab = ylab, as.table = TRUE, ...))
            }
            else {
                ## plotting of both points and lines not compatible with group
                if (nt > 1) {
                    key <- list(space = space, corner = corner,
                                lines = list(col = col, lty = lty, lwd = lwd,  pch = pch),
                                type = "b", text = list(lt), border = TRUE,
                                columns = columns)
                }
                else {
                    key <- NULL
                    pch <- rep(pch[1], nlevels(treatment))
                    col <- rep(col[1], nlevels(treatment))
                }
                ylim <- range(c(x$xscores, x$yscores))
                print(xyplot(x$xscores ~ x$time | subject, group=treatment,
                             type=c("l"), col = col, lty = lty, lwd = lwd,
                             main = ifelse(missing(main),
                                           'Fitted Values Corresponding to Different Subjects', main),
                             xlab = xlab,  ylab = ylab, as.table = TRUE, key = key,
                             ylim = ylim + diff(ylim)/20*c(-1, 1), ...))
                if (scores) {
                    layout <- trellis.currentLayout()
                    nc <- ncol(layout)
                    for (i in seq_len(ns)) {
                        trellis.focus("panel", column = (i - 1) %% nc + 1,
                                      row = (i - 1) %/% nc + 1)
                        id1 <- subject == ls[i]
                        id2 <- unclass(as.factor(treatment))[id1]
                        panel.points(x$time[id1], x$yscores[id1], col = col[id2],
                                     pch = pch[id2], cex = 0.6)
                        trellis.unfocus()
                    }
                }
            }
        }
        else {
            if (signature) {
                for (i in seq_len(ns)) {
                    plot(freq, x$ycoef[,i],xlab=xlab,ylab=ylab,type='l',
                         col = col[1], lty = lty, lwd = lwd, xaxt = "n", ...)
                    axis(1, xaxis$at, xaxis$freq.nm)
                    title(paste(ifelse(missing(main), 'Signature', main),
                                paste(" ", ls[i])[!one]))
                }
            }
            else {
                yscores <- tapply(x$yscores, list(x$time, treatment, subject), I)
                xscores <- tapply(x$xscores, list(x$time, treatment, subject), I)
                for (i in seq_len(ns)) {
                    matplot(names(yscores[, 1, i]), yscores[, , i],
                            col = col, type = ifelse(scores, "p", "n"), pch = pch,
                            xlab = xlab,  ylab = ylab, cex = 0.6, ...)
                    matlines(names(xscores[, 1, i]), xscores[, , i], col = col,
                             type = "l", lty = lty, lwd = lwd)
                    title(paste(ifelse(missing(main), 'Fitted Values', main),
                                paste(" ", ls[i])[!one]))
                    if (nt > 1) {
                        legend(x = legend.x, y=NULL, legend = lt, col = col,
                               lty = lty, lwd = lwd, pch = pch, merge = TRUE)
                    }
                }
            }
        }
    }
    else {
        if (signature) {
            plot(rowMeans(x$ycoef),xlab=xlab, ylab = ylab, type='l',
                 col = col[1], lty = lty, lwd = lwd, xaxt = "n", ...)
            axis(1, xaxis$at, xaxis$freq.nm)
            title(ifelse(missing(main), 'Mean Signature', main))
        }
        else {
            xscores <- tapply(x$xscores, list(x$time, treatment), mean)
            matplot(rownames(xscores), xscores, col = col, type = "l", lty = lty, lwd = lwd,
                    xlab = xlab,  ylab = ylab, ...)
            title(ifelse(missing(main), 'Mean Fitted Values', main))
            if (nt > 1) {
                legend(x = legend.x,y=NULL, legend = lt, col = col, lty = lty, lwd = lwd,
                       pch = NULL, merge = TRUE)
            }
       }
    }
    invisible()
}
