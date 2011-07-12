print.gslcca <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

    if (x$global.roots != 0)
        cat("Data pre-smoothed using", x$global.roots, "roots\n")

    ns <- ncol(x$ycoef)
    if (ns != 1) {
        cat("GSLCCA based on", ns, "subjects\n")
        if (x$subject.roots != 0) {
            cat("\nData smoothed separately for each subject using", x$subject.roots,
                "roots\n")
        }
    }

    cat("\nNonlinear parameters:\n")
    print(format(x$nonlinear.parameters, digits = max(4, digits + 1)))
}
