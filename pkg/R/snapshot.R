snapshot <- function(x, y, 
                     plot = c("signatures", "boxplot", "mean", "pvalue"),
                     names = NULL, normalise = FALSE, 
                     breaks = c(4, 8, 13, 30),
                     col = NULL, lty = 1, lwd = 1.5, ...){
    if (nrow(x) != nrow(y)){
        nr <- min(nrow(x), nrow(y))
        x <- x[1:nr,]
        y <- y[1:nr,]
    }
    if (normalise){
        norm <- function(x) x/sqrt(sum(x^2))
        x <- apply(x, 2, norm)
        y <- apply(y, 2, norm)
    }
    
    freq <- as.numeric(rownames(x))
    ns <- ncol(x)
    
    if (missing(col)) {
        if (exists("rainbow_hcl", mode= "function"))
            col <- rainbow_hcl(ns, c = 100, l = 60) 
        else
            col <- seq_len(ns)
    }
    
    signatures <- length(agrep("signatures", plot))
    boxplot <- length(agrep("boxplot", plot))
    mean <- length(agrep("mean", plot))
    manhattan <- length(agrep("pvalue", plot))
    
    ### signatures
    if (signatures) {
        matplot(freq, x, col=col, type='l', xlab= "Frequency (Hz)", 
                ylab= "Coefficient", lty = lty, lwd = lwd, ...)
        title(paste("Signatures for", names[1]))
        matplot(freq, y, col=col, type='l', xlab= "Frequency (Hz)", 
                ylab= "Coefficient", lty = lty, lwd = lwd, ...)
        title(paste("Signatures for", names[2]))
    }
    
    ### boxplot
    if (boxplot) {
        boxplot(value ~ X2, data = melt(t(x)), ylab = "Loading", 
                xlab = "Frequency (Hz)", main = names[1])
        boxplot(value ~ X2, data = melt(t(y)), ylab = "Loading", 
                xlab = "Frequency (Hz)", main = names[2])
    }
    
    ### plot mean signatures
    if (mean){
        y1 <- rowMeans(x)
        y2 <- rowMeans(y)
        ylim <- range(y1, y2)
        plot(seq_along(y1), y1, type = "b", ylab = "Signature", 
            xlab = "Frequency (Hz)", col = "red", ylim = ylim)
        abline(h = 0)
        lines(seq_along(y1), y2, col = "blue", type = "b")
        legend("topright", y=NULL, legend = names, col = c("red", "blue"),
               lty = lty, lwd = lwd, pch = NULL, merge = TRUE) 
    }

    ## do t-tests for each row of x vs each row of y
    grp <- as.factor(rep(1:2, c(ncol(x), ncol(y))))
    pvalue <- rowttests(as.matrix(cbind(x, y)), grp)$p.value

    ## plot Manhattan plot
    if (manhattan) {
        xx <- seq_along(y1)
        yy <- p.adjust(pvalue, method = "fdr")
        yy <- -log(yy, 10)
        ylim <- range(c(yy, -log(0.05, 10)))
        plot(xx, yy, type="n", xlab="Frequency", ylab=expression(-log[10](p)),
             ylim = ylim)
        col <- c(hex(polarLUV(H = 0, C = 0, L = 60)),
                 rainbow_hcl(4, c = 100, l = 60))
        breaks1 <- c(1, breaks)
        breaks2 <- c(breaks, ceiling(max(freq)))
        bands <- mapply(seq, breaks1, breaks2)
        for (b in 1:5){
            x <- c(xx[bands[[b]][1]], xx[bands[[b]]],
                   xx[rev(bands[[b]])[1]])
            y <- c(0, yy[bands[[b]]], 0)
            polygon(x, y, col = col[b], border = col[b])
        }
        points(xx, yy, pch = 4)
        lines(xx, yy, pch = 4)
        abline(h = -log(0.05, 10), lty = 2) #significance line
        legend("topright", c("delta", "theta", "alpha", "beta", "gamma"), 
             col = col[1:5], pch = 15)
        }
    invisible(pvalue)
}