summary.ESLCCA <- function(object, ...){
    res <- list(call = object$call,
                global.roots = object$global.roots,
                subject.roots = object$subject.roots,
                nonlinear.parameters = object$nonlinear.parameters,
                mean.signature = rowMeans(as.matrix(object$ycoef)),
                mean.fitted = tapply(object$xscores, list(object$time, treatment), mean))
    class(res) <- "summary.ESLCCA"
    res
}
