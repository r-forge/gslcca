varySmooth <- function(x, subject.smooth = 1:10, ...){
    if(!inherits(x, "gslcca")) stop("'x' must be an \"gslcca\" object")

    ## fit models varying subject.smooth
    res <- list()
    for (r in subject.smooth) {
        if (r == x$subject.roots) res[[as.character(r)]] <- x
        else {
            cat("Re-running gslcca with subject.smooth =", r, "\n")
            res[[as.character(r)]] <- update(x, subject.smooth = r)
        }
    }

    attributes(res) <- list(subject.smooth = subject.smooth, class = "varySmooth")
    res
}
