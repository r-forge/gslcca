print.gslcca <- function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

    if (x$global.smooth != 0)
        cat("Data pre-smoothed using", x$global.smooth, "roots\n")

    ns <- ncol(x$ycoef)
    if (ns != 1) {
        cat("GSLCCA based on", ns, "subjects\n")
    }
    if (x$subject.smooth != 0) {
        cat("\nData smoothed at subject level using", x$subject.smooth,
            "roots\n")
    }

    cat("\nNonlinear parameters:\n")
    print(format(x$nonlinear.parameters, digits = max(4, digits + 1)),
          quote = FALSE)
    
}
