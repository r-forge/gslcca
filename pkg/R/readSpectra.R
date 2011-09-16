readSpectra <- function(file, 
                        info = 1:3,
                        treatment = c("Control", "Low Dose", 
                                      "Middle Dose", "High Dose"),
                        resolution = 1,
                        end = 42900,
                        nfreq = 36,
                        ...) {
    ## read in data
    dat <- read.delim(file, check.names = FALSE, ...)
    
    ## split out spectra
    spectra <- as.matrix(dat[ , -info])
    dat <- dat[,info] 
    dat$spectra <- spectra
    
    ## subset up to end time
    nm <- tolower(colnames(dat))
    dat <- dat[dat[,nm == "time"] <= end,]
    
    ## fix treatment levels
    if ("treatment" %in% nm) { 
        lev <- levels(dat[,nm == "treatment"])
        if (!setequal(lev, treatment)) 
            stop("Valid levels of treatment factor specified as",
                 paste(treatment, sep = " "),
                "but unique values found to be",
                 paste(lev, sep = " "))
        dat[,nm == "treatment"] <- factor(dat[,nm == "treatment"], 
                                          levels  = treatment)
    }
    
    ## average over frequencies at desired resolution
    freq <- as.numeric(colnames(dat$spectra)) %/% resolution
    if (any(is.na(freq)))
        stop("Column names of spectra matrix must be numeric")
    if (any(duplicated(freq))){
        d <- as.numeric(names(sort(-table(tabulate(freq))))[1]) 
        ## average every d frequencies
        bands <- rep(1:(length(freq)/d), each = d)
        nm <- colnames(dat$spectra)
        dat$spectra <- t(rowsum(t(dat$spectra), bands, 
                                reorder = FALSE, na.rm = TRUE))
        colnames(dat$spectra) <- nm[!as.logical(seq_along(freq) %% d)]
    }
    dat$spectra <- dat$spectra[, 1:(nfreq%/%resolution)]
    dat
}