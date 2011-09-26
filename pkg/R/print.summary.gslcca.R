print.summary.gslcca <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), sep = "", fill = TRUE)

    if (x$global.smooth != 0)
        cat("Data pre-smoothed using", x$global.smooth, "roots\n\n")

    ns <- NCOL(x$ycoef)
    if (ns != 1) {
        cat("gslcca based on", ns, "subjects\n")
    }
    if (x$subject.smooth != 0) {
        cat("\nData smoothed at subject level using", x$subject.smooth,
            "roots\n")
    }
    cat("\nPercent variance explained by SVD approximation:\n")
    print(format(x$pct.explained, digits = digits), quote = FALSE)

    cat("\nNonlinear parameters:\n")
    print(format(x$nonlinear.parameters, digits = digits), quote = FALSE)

    cat("\nCorrelation at final iteration:\n")
    print(format(x$cor, digits = digits), quote = FALSE)

    if (any(x$opt.conv)) {
        not.conv  <- names(x$opt.conv)[as.logical(x$opt.conv)]
         cat("\nThe algorithm did not converge",
             "for the following subjects:\n", not.conv)
         cat("\nThe maximum number of iterations was ", max(x$opt.iter), ".\n",
             sep = "")
    }
}
