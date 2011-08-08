plot.varySmooth <- function(x, type = "opt", series = x[[1]]$treatment, subject = levels(x[[1]]$subject),
                            ask = dev.interactive(), main = NULL, xlab = NULL, ylab = NULL,
                            col = NULL, lty = 1, lwd = 1.5, pch = NULL,
                            space = "bottom", corner = NULL, columns = 2, ...){
    if(!inherits(x, "varySmooth")) stop("'x' must be an object of class \"varySmooth\"")

    subject.smooth <- attr(x, "subject.smooth")
    nr <- length(subject.smooth)

    subj <- x[[1]]$subject
    if (is.null(subj)) subj <- factor(1)
    ns <- nlevels(subj)

    treatment <- series
    nt  <- nlevels(treatment)
    if (is.null(treatment)) {
        ## create factor to id series
        time <- as.numeric(paste(unclass(subj), x[[1]]$time, sep = ""))
        ord <- order(time)
        reps <- table(time)
        treatment <- factor(sequence(reps)[order(ord)])
        nt <- 1
    }

    if (length(agrep("opt", type))) type <- "opt"
    else if (length(agrep("fitted", type))) type <- "fitted"
    else if (length(agrep("scores", type))) type <- "scores"
    else if (length(agrep("signature", type))) type <- "signature"
    else stop("'type' not one of recognised options: \"opt\", \"fitted\", \"scores\", or \"signature\"")

    ## sort colours and lines
    if (missing(xlab)) xlab <- switch(type,
                                      signature = "Frequency (Hz)",
                                      fitted = "Time",
                                      opt = "Roots")
    if (missing(ylab)) ylab <- switch(type,
                                      signature = "Coefficient",
                                      fitted = "Score",
                                      opt = expression(Log(1 - R^2)))
    n <- ifelse(type == "signature", ns, nt)
    if (missing(col)) {
        if (exists("rainbow_hcl", mode= "function"))
            col <- c(hex(polarLUV(H = 0, C = 0, L = 60)), rainbow_hcl(n - 1, c = 100, l = 60)) # assumes ref
        else
            col <- seq_len(n)
    }
    if (missing(pch)) pch <- seq_len(n)


    if (type == "opt") {
        opt <- lapply(x, function(x) sapply(x$opt, "[[", "value"))
        roots <- rep(subject.smooth, each = length(opt[[1]]))
        if (!is.null(x[[1]]$subject)){
            subj <- gl(length(opt[[1]]), 1, length(roots), labels = levels(x[[1]]$subject))
            f <- unlist(opt) ~ roots | subj
        }
        else f <- unlist(opt) ~ roots
        miny <- min(unlist(opt))
        print(xyplot(f, type="b", main = ifelse(missing(main), "Optimized Value", main), xlab= "Roots",
                     ylab = expression(Log(1 - R^2)), as.table = TRUE, ylim = c(miny + miny/20, 0), ...))
    }

    if (type == "fitted" | type == "scores") {
        ##need to do separately for each subject
        if (is.null(subject)) subject <- levels(subj)
        else if (any(!subject %in% levels(subj)))
            stop("the following values of 'subject' are not valid: ", setdiff(subject, levels(subj)))
        if (length(subject) == 1) ask <- FALSE
        if (ask) {
            oask <- devAskNewPage(NULL)
            on.exit(devAskNewPage(oask))
        }

        for (s in subject) {
            main0 <- paste("Fitted Values", ifelse(is.null(x[[1]]$subject), "", paste("Subject", s)))
            xscores <- sapply(x, function(x) subset(x$xscores, subj == s))
            if (type == "scores") yscores <- c(sapply(x, function(x) subset(x$yscores, subj == s)))
            time <- c(sapply(x, function(x) subset(x$time, subj == s)))
            trt <- rep(subset(treatment, subj == s), nr)
            roots <- gl(nr, nrow(xscores), labels = subject.smooth)
            if (nt > 1) {
                key <- list(space = space, corner = corner, lines = list(col = col, lty = lty, lwd = lwd,  pch = pch),
                            type = "b", text = list(levels(treatment)), border = TRUE, columns = columns)
            }
            else {
                key <- NULL
                pch <- rep(pch[1], nlevels(treatment))
                col <- rep(col[1], nlevels(treatment))
            }
            if (ask) devAskNewPage(TRUE)
            if (type == "scores") ylim <- range(c(xscores, yscores))
            else ylim <- range(xscores)
            print(xyplot(c(xscores) ~ time | roots, group = trt, type = "l", col = col,
                         lty = lty, lwd = lwd, main = ifelse(missing(main), main0, main), ylim = ylim + diff(ylim)/20*c(-1, 1),
                         xlab = "Time",  ylab = "Score", as.table = TRUE, key = key, ...))
            devAskNewPage(FALSE)
            if (type == "scores") {
                layout <- trellis.currentLayout()
                nc <- ncol(layout)
                for (i in seq_len(nr)) {
                    trellis.focus("panel", column = (i - 1) %% nc + 1, row = (i - 1) %/% nc + 1)
                    id1 <- roots == subject.smooth[i]
                    id2 <- unclass(as.factor(trt[id1]))
                    panel.points(time[id1], yscores[id1], col = col[id2],
                                 pch = pch[id2], cex = 0.6)
                    trellis.unfocus()
                }
            }
        }
    }

    if (type == "signature") {
        ycoef <- lapply(x, "[[", "ycoef")
        nf <- nrow(ycoef[[1]])
        freq <- rep(seq_len(nf), ns * nr)
        roots <- gl(nr, ns * nf, labels = subject.smooth)
        if (!is.null(x[[1]]$subject))
            subj <- gl(ns, nf, length(roots), labels = levels(x[[1]]$subject))
        else subj <- NULL
        if (ns > 1) {
            key <- list(space = space, corner = corner, lines = list(col = col, lty = lty, lwd = lwd),
                        type = "l", text = list(levels(subj)), border = TRUE, columns = columns)
        }
        else key <- NULL
        xyplot(unlist(ycoef) ~ freq | roots, group=subj, type="l",
               main = ifelse(missing(main), 'Signature by Number of Roots', main),
               col = col, lty = lty, lwd = lwd, xlab = xlab,  ylab = ylab, as.table = TRUE, key = key, ...)
    }
}
