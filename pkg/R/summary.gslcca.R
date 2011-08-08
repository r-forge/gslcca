summary.gslcca <- function(object, ...){
    opt.value <- sapply(object$opt, "[[", "value")
    opt.iter <- sapply(object$opt, function(x) x$counts[1])
    opt.conv <- sapply(object$opt, "[[", "convergence")
    names(opt.value) <- names(opt.iter) <- names(opt.conv) <- colnames(object$xcoef)
    res <- list(call = object$call,
                global.roots = object$global.roots,
                subject.roots = object$subject.roots,
                nonlinear.parameters = object$nonlinear.parameters,
                cor = object$cor,
                opt.value = opt.value,
                opt.iter = opt.iter,
                opt.conv = opt.conv)
    class(res) <- "summary.gslcca"
    res
}
