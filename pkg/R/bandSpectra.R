bandSpectra <- function(spectra, breaks = NULL, labels = NULL, ...){
    ## combine frequencies
    freq <- as.numeric(colnames(spectra))
    breaks <- c(floor(min(freq - 1, 0)), breaks, ceiling(max(freq) + 1))
    bands <- cut(freq, breaks = breaks, labels = labels, ...)   
    
    total <- t(rowsum(t(spectra), bands, reorder = FALSE, na.rm = TRUE))
    count <- t(rowsum(t(1 - is.na(spectra)), bands, reorder = FALSE, na.rm = TRUE))
    total/count
}