bandSpectra <- function(spectra, breaks = NULL, labels = NULL, ...){
    ## convert breaks to frequency bands
    freq <- as.numeric(colnames(spectra))
    breaks <- c(floor(min(freq - 1, 0)), breaks, ceiling(max(freq) + 1))
    bands <- cut(freq, breaks = breaks, labels = labels, ...)   
    
    ## total power in each band
    t(rowsum(t(spectra), bands, reorder = FALSE, na.rm = TRUE))
}