print.summary.gslcca <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), sep = "", fill = TRUE)

    if (x$global.roots != 0)
        cat("Data pre-smoothed using", x$global.roots, "roots\n\n")

    ns <- NCOL(x$ycoef)
    if (ns != 1) {
        cat("gslcca based on", ns, "subjects\n")
        if (x$subject.roots != 0) {
            cat("\nData smoothed separately for each subject using", x$subject.roots,
                "roots\n")
        }
    }

    cat("\nNonlinear parameters:\n")
    print(format(x$nonlinear.parameters, digits = digits))

    cat(ifelse(ns != 1, "\nMean signature:\n", "\nSignature:\n"))
    print.default(format(x$mean.signature, digits = digits), quote = FALSE)

    cat(ifelse(ns != 1, "\nMean fitted values:\n", "\nFitted values:\n"))
    print(format(x$mean.fitted, digits = digits), quote = FALSE)
}
