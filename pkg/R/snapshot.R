snapshot <- function(x, y, 
                     p.adjust.method = "fdr",
                     normalise = TRUE, 
                     ...){  
    nr <- min(nrow(x), nrow(y))    
    if (nrow(x) != nrow(y)){
        x <- x[1:nr,]
        y <- y[1:nr,]
    }
    if (normalise){
        norm <- function(x) x/sqrt(sum(x^2))
        x <- apply(x, 2, norm)
        y <- apply(y, 2, norm)
    }   

    ## do t-tests for each row of x vs each row of y
    grp <- as.factor(rep(1:2, c(ncol(x), ncol(y))))
    pvalue <- p.adjust(rowttests(as.matrix(cbind(x, y)), grp, ...)$p.value, 
                       method = p.adjust.method)

    res <- list(call = match.call(), x = x, y = y, pvalue = pvalue)
    class(res) <- "snapshot"
    res
}