print.summary.gslcca <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), sep = "", fill = TRUE)

    if (x$global.roots != 0)
        cat("Data pre-smoothed using", x$global.roots, "roots\n\n")

    ns <- NCOL(x$ycoef)
    if (ns != 1) {
        cat("gslcca based on", ns, "subjects\n")
    }
    if (x$subject.roots != 0) {
        cat("\nData smoothed at subject level using", x$subject.roots,
            "roots\n")
    }
    cat("\nPercent variance explained by SVD approximation:\n")
    print(format(x$pct.explained, digits = digits))

    cat("\nNonlinear parameters:\n")
    print(format(x$nonlinear.parameters, digits = digits))

    cat("\nCorrelation at final iteration:\n")
    print(format(x$cor, digits = digits), quote = FALSE)

    cat("\nNumber of iterations:\n")
    print(format(x$opt.iter, digits = digits), quote = FALSE)

    cat("\nConvergence indicator (0 = converged, 1 = not converged):\n")
    print(format(x$opt.conv, digits = digits), quote = FALSE)
}
