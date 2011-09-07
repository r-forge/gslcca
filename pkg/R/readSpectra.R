readSpectra <- function(file, 
                        info = 1:3,
                        treatment = c("Control", "Low Dose", 
                                      "Middle Dose", "High Dose"),
                        resolution = 1,
                        ...) {
    ## read in data
    dat <- read.delim(file, check.names = FALSE, ...)
    
    ## split out spectra
    spectra <- as.matrix(dat[ , -info])
    dat <- dat[,info]    
    dat$spectra <- spectra
    
    ## fix treatment levels
    nm <- tolower(colnames(dat))
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
    
    ## keep only frequencies at desired resolution
    freq <- as.numeric(colnames(dat$spectra)) %/% resolution
    if (any(is.na(freq)))
        stop("Column names of spectra matrix must be numeric")
    if (any(duplicated(freq))){
        d <- as.numeric(names(sort(-table(tabulate(freq))))[1]) 
        ## want every d'th freq
        dat$spectra <- dat$spectra[ , !as.logical(seq_along(freq) %% d)]
    }
    dat
}