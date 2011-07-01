varySmooth <- function(x, subject.smooth = NULL, ...){
    if(!inherits(x, "ESLCCA")) stop("'x' must be an \"ESLCCA\" object")

    if (is.null(subject.smooth)) subject.smooth <- seq_len(ncol(x$ycoef)/2)

    ## fit models varying subject.smooth
    res <- list()
    for (r in subject.smooth) {
        if (r == x$subject.roots) res[[as.character(r)]] <- x
        else {
            cat("Re-running ESLCCA with subject.smooth =", r, "\n")
            res[[as.character(r)]] <- update(x, subject.smooth = r)
        }
    }

    attributes(res) <- list(subject.smooth = subject.smooth, class = "varySmooth")
    res
}
