snapshot.mod <- function(x, y, compact = TRUE, p.adjust.method = "fdr",
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

   ### plot mean signatures
    plot.mean <- function(x, y, freq, col, pch.col, pch, legend, lines = TRUE, 
                          add = FALSE, ...){
        y1 <- rowMeans(x)
        y2 <- rowMeans(y)
        ylim <- range(y1, y2)
        if (!add) {
            plot(freq, y1, type = "n", ylab = "", 
                 xlab = "", ylim = ylim, ...)
            abline(h = 0)
        }
        if (lines){
            if (!is.list(col)) col <- list(col, col)
            if (!is.list(pch)) pch <- list(pch, pch)
            if (!is.list(pch.col)) pch.col <- list(pch.col, pch.col)
            lines(freq, y1, col = col[[1]], lwd = 2, ...) 
            lines(freq, y2, col = col[[2]], lwd = 2, ...)    
            points(freq, y1, col = pch.col[[1]], pch = pch[[1]], ...) 
            points(freq, y2, col = pch.col[[2]], pch = pch[[2]], ...)
        }
        if (!missing(legend)) {
            legend("topright", legend, col = unlist(col), 
                   pch = unlist(pch))
        }
    }

    ## do t-tests for each row of x vs each row of y
    grp <- as.factor(rep(1:2, c(ncol(x), ncol(y))))
    pvalue <- p.adjust(rowttests(as.matrix(cbind(x, y)), grp)$p.value, 
                       method = p.adjust.method)

    col <- rainbow_hcl(length(breaks) + 1, c = 35, l = 85)
    col <- col[c(seq(1, length(col), 2), seq(2, length(col), 2))]

        cutoff <- cut(pvalue, c(1, 0.05, 0.01, 0.001, 0))
        cols <- c(heat_hcl(4, h= c(0, 100), c = c(200, 30), 
                           l = c(40, 90), power = c(1/5, 3))[1:3], 
                  hex(polarLUV(H = 0, C = 0, L = 70)))
        plot.mean(x, y, freq, lines = FALSE, ...)    
        usr <- par("usr")
        breaks1 <- c(-1, breaks)
        breaks2 <- c(breaks, ceiling(usr[2]) + 1)
        poly <- mapply(seq, breaks1, breaks2)              
        for (b in seq_along(poly)){
            xx <- c(poly[[b]][1], poly[[b]],
                   rev(poly[[b]])[1])
            r <- length(xx) - 2
                polygon(xx, c(usr[3], rep(usr[4], r), usr[3]), 
                    col = col[b], border = col[b])
        }
        axis(3, at = breaks, labels = FALSE, tick = FALSE)
        mtext(bands, side = 3, 
              at = diff(c(usr[1], breaks, usr[2]))/2+ c(0, breaks), 
              line = 0.5, cex = 1.8)
        plot.mean(x, y, freq, pch.col = cols[cutoff], pch = list(15, 17), 
                  col = cols[4], add = TRUE, ...)                
        legend("topright", bty = "n", 
               c(names, "p < 5%", "p < 1%", "p < 0.1%"), 
               col = c(cols[4], cols[4], cols[3:1]), 
               pch = c(15, 17, 15, 15, 15), cex = 1.3)            
}