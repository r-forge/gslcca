snapshot <- function(x, y, 
                     plot = c("signatures", "mean", "pvalue", "compact"),
                     p.adjust.method = "fdr",
                     names = NULL, normalise = FALSE, 
                     breaks = c(4, 8, 13, 30),
                     bands  =  c("delta", "theta", "alpha", "beta",
                                 "gamma"),
                     col = NULL, lty = 1, lwd = 1.5, ...){  
    freq <- as.numeric(rownames(x))
    ns <- ncol(x)
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
    
    if (missing(col)) {
        if (exists("rainbow_hcl", mode= "function"))
            col <- rainbow_hcl(ns, c = 100, l = 60) 
        else
            col <- seq_len(ns)
    }
    
    signatures <- length(agrep("signatures", plot))
    mean <- length(agrep("mean", plot))
    manhattan <- length(agrep("pvalue", plot))
    compact <- length(agrep("compact", plot)) 
    
    M <- list(x, y)
    
    ### signatures
    if (signatures) {
        for (i in 1:2) {
            matplot(freq, M[[i]], col=col, type='l', 
                    xlab= "Frequency (Hz)", ylab= "Coefficient", 
                    lty = lty, lwd = lwd, ...)
            title(paste("Signatures for", names[i]))
        }
    }
    
    ### plot mean signatures
    plot.mean <- function(x, y, freq, col, pch.col, pch, legend, lines = TRUE, 
                          add = FALSE, ...){
        y1 <- rowMeans(x)
        y2 <- rowMeans(y)
        ylim <- range(y1, y2)
        if (!add) {
            plot(freq, y1, type = "n", ylab = "Signature", 
                 xlab = "Frequency (Hz)", ylim = ylim, ...)
            abline(h = 0)
        }
        if (lines){
            if (!is.list(col)) col <- list(col, col)
            if (!is.list(pch)) pch <- list(pch, pch)
            if (!is.list(pch.col)) pch.col <- list(pch.col, pch.col)
            lines(freq, y1, col = col[[1]]) 
            lines(freq, y2, col = col[[2]])    
            points(freq, y1, col = pch.col[[1]], pch = pch[[1]]) 
            points(freq, y2, col = pch.col[[2]], pch = pch[[2]])
        }
        if (!missing(legend)) {
            legend("topright", legend, col = unlist(col), 
                   pch = unlist(pch))
        }
    }
            
    if (mean){
        plot.mean(x, y, freq, col = list("red", "blue"), pch = 1, 
                  pch.col = list("red", "blue"), legend = names)
    }

    ## do t-tests for each row of x vs each row of y
    grp <- as.factor(rep(1:2, c(ncol(x), ncol(y))))
    pvalue <- p.adjust(rowttests(as.matrix(cbind(x, y)), grp)$p.value, 
                       method = p.adjust.method)

    col <- rainbow_hcl(length(breaks) + 1, c = 60, l = 75)
    col <- col[c(seq(1, length(col), 2), seq(2, length(col), 2))]
      
    if (compact) {
        cutoff <- cut(pvalue, c(1, 0.05, 0.01, 0.001, 0))
        cols <- heat_hcl(4, c. = c(100, 30), l = c(50, 90))
        plot.mean(x, y, freq, lines = FALSE, ...)    
        usr <- par("usr")
        breaks1 <- c(-1, breaks)
        breaks2 <- c(breaks, ceiling(usr[2]) + 1)
        poly <- mapply(seq, breaks1, breaks2)              
        for (b in seq_along(poly)){
            xx <- c(poly[[b]][1], poly[[b]],
                   rev(poly[[b]])[1])
            r <- length(xx) - 2
            yy <- list (c(-0.002, rep(0.002, r), -0.002),
                       c(usr[4] -0.004, rep(usr[4], r), usr[4] -0.004),
                       c(usr[3], rep(usr[3] + 0.004, r), usr[3]))
            for (i in 1:3)
                polygon(xx, yy[[i]], col = col[b], border = col[b])
        }
        axis(3, at = breaks, labels = FALSE)
        mtext(bands, side = 3, 
              at = diff(c(usr[1], breaks, usr[2]))/2+ c(0, breaks))
        plot.mean(x, y, freq, pch.col = cols[cutoff], pch = list(15, 17), 
                  col = cols[4], add = TRUE)                
        legend("topright", bty = "n", 
               c(names, "p < 5%", "p < 1%", "p < 0.1%"), 
               col = c(cols[4], cols[4], cols[3:1]), 
               pch = c(15, 17, 15, 15, 15))
    }

    ## plot Manhattan plot
    if (manhattan) {
        xx <- seq_len(nr)
        yy <- -log(pvalue, 10)
        ylim <- range(c(yy, -log(0.05, 10)))
        plot(xx, yy, type="n", xlab="Frequency", 
             ylab=expression(-log[10](p)), ylim = ylim)
        breaks1 <- c(1, breaks)
        breaks2 <- c(breaks, ceiling(max(freq)))
        poly <- mapply(seq, breaks1, breaks2)
        for (b in 1:5){
            x <- c(xx[poly[[b]][1]], xx[poly[[b]]],
                   xx[rev(poly[[b]])[1]])
            y <- c(0, yy[poly[[b]]], 0)
            polygon(x, y, col = col[b], border = col[b])
        }
        points(xx, yy, pch = 4)
        lines(xx, yy, pch = 4)
        abline(h = -log(0.05, 10), lty = 2) #significance line
        legend("topright", bands, col = col, pch = 15)
    }
    
    invisible(pvalue)
}