summary.ESLCCA <- function(object, ...){
    list(call = x$call,
         global.roots = x$global.roots,
         subject.roots = x$subject.roots,
         nonlinear.parameters = x$nonlinear.parameters,
         mean.signature = colMeans(x$ycoef),
         mean.fitted = cast(x$scores, time ~ treatment, value = "xscores", fun.aggregate = mean))
}
