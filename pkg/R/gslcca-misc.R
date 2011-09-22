signatures <- function(object, ...){
    if(!inherits(object, "gslcca")) 
        stop("'object' must be an object of class \"gslcca\"")
    object$ycoef
}

fitted.gslcca <- function(object, ...){
    object$xscores
}