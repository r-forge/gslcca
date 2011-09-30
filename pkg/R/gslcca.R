gslcca <- function (Y, # matrix of power spectra
                    formula = "Double Exponential",
                    time, # time points
                    subject = NULL,
                    global = FALSE, ##fit one model for all subjects?
                    treatment = NULL,
                    ref = 1,
                    separate = FALSE,# use "proportional" instead c.f. hazard?
                    partial = ~1,
                    data = NULL,
                    subset = NULL,
                    global.smooth = FALSE,
                    subject.smooth = TRUE,
                    pct.explained = 0.96,
                    start = NULL,
                    method = "L-BFGS-B",
                    lower = 2,
                    upper = 15,
                    ...) { # allow arguments to be passed to optim
    ## Formula must be a function of "time", possibly other var and parameters
    if (is.character(formula)){
        models <- c("Double Exponential", "Critical Exponential")
        formula <- models[agrep(formula, models)]
        if (!length(formula))
            stop("\"formula\" not recognised, only \"Double Exponential\"",
                 "or \"Critical Exponential\" \n",
                "can be specified as character strings")
        if (is.null(start))
            start <- switch(formula, #replicate later for multiple treatments
                            #K1 > K2 so increases from ref
                            "Double Exponential" = list(K1 = 9, K2 = 8.5),
                            "Critical Exponential" = list(K1 = 8.5))
        formula <- switch(formula,
                          "Double Exponential" = ~ exp(-time/exp(K1))-exp(-time/exp(K2)),
                          "Critical Exponential" = ~ time*exp(-time/exp(K1)))
    }

    vars <- all.vars(formula)
    if (!"time" %in% vars) stop("'formula' must be a function of 'time'")

    ## anything without starting values assumed to be a variable
    vars <- setdiff(c(vars, all.vars(partial)), c("time", names(start)))
    dummy <- reformulate(c("0", vars)) #don't need int for mf
    formula <- as.expression(formula[[length(formula)]])

    ## Collate data to ensure equal length and deal with NAs
    ## get model frame inc. Y, time, subject & treatment if specified, omit NAs
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("Y", "subject", "treatment", "time", "data", "subset"),
               names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- dummy
    mf$na.action <- na.omit
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## split up again for convenience
    Y <- mf$`(Y)`
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    nr <- nrow(Y)
    freq.nm <- colnames(Y)
    f <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", freq.nm)))
    if (any(is.na(f))) f <- seq_along(f)

    if (!is.null(mf$`(subject)`)) {
        subject <- as.factor(mf$`(subject)`)
        ind <- split(seq_len(nr), subject)
    }
    else {
        subject <- NULL
        ind <- list(seq_len(nr))
    }
    nind <- length(ind)

    if (!is.null(mf$`(treatment)`)) {
        treatment <- as.factor(mf$`(treatment)`)
        if (!is.numeric(ref)) ref <- which(levels(treatment) == ref)
        ntrt <- nlevels(treatment) - 1
    }
    else ntrt <- 1

    if (is.null(mf$`(time)`))
        stop("'time' must be specified and have length equal to", nr, "\n")
    names(mf)[match("(time)", names(mf))] <- "time"

    ## take out extras from mf
    mf <- mf[, !(names(mf) %in% c("(Y)", "(subject)", "(treatment)")),
             drop = FALSE]

    ## Pre-smoothing of the complete data set for artefacts removal
    if (global.smooth) {
        if (global.smooth == TRUE) {
            Ytemp= scale(Y,scale=FALSE)
            eig=eigen(crossprod(Ytemp), only.values = TRUE)$values
            nroots=seq_along(eig)
            global.smooth=max(nroots[(eig/eig[1])>0.001 &
                cumsum(eig)/sum(eig)<pct.explained]) + 1
        }
        SVD = svd( Y, nu=global.smooth, nv=global.smooth )
        Y = SVD$u%*% (SVD$d[1:global.smooth] * t(SVD$v))
    }

    cum.pct.explained <- matrix(, nrow = ncol(Y), ncol = nind)
    r <- numeric(nind)
    for (i in seq_len(nind)) {
        Ytemp = scale(Y[ind[[i]],],scale=FALSE)
        eig = eigen(crossprod(Ytemp), only.values = TRUE)$values
        cum.pct.explained[,i] <- cumsum(eig)/sum(eig)
        r[i] = sum((eig/eig[1])>0.001 & cum.pct.explained[,i] < pct.explained)
    }
    if (isTRUE(subject.smooth)) {
        ## Set number of roots to satisfy pct.explained
        subject.smooth <- max(r)+1
    }
    else if(!is.numeric(subject.smooth)) subject.smooth <- ncol(Y)
    cum.pct.explained <- cum.pct.explained[subject.smooth,]

    reps <- ifelse(separate, ntrt, 1)
    if (is.null(start)) {
        stop('No initial values have been specified for the nonlinear',
             ' parameters.')
    }
    else {
        start <- mapply(function(start, label, reps) {
            if (length(start) == 1) rep(start, reps)
            else if (length(start) == reps) start
            else stop("there should be 1 or", reps,
                      "starting values for parameter", label, "\n")
        }, start = start, label = names(start), MoreArgs = list(reps = reps),
                        SIMPLIFY = FALSE)
    }
    nPar <- length(start)
    group <- gl(nPar, reps, labels = names(start))

    ## evaluate RHS for first row of data
    par <- as.list(rep(NA, length(start)))
    names(par) <- names(start)
    dat <- do.call("cbind", list(mf[1,, drop = FALSE], par))
    val <- eval(formula, dat)
    mt <- terms(partial)
    CCA.roots = 1

    ny <- length(f)*nind
    if (nind > 1 & global == TRUE){
        nind <- 1
        ind <- list(seq_len(nr))
    }
    ycoef <- matrix(numeric(ny), ncol=nind)    
    xcoef <- matrix(nrow=length(val)*ntrt, ncol=nind)
    nonlin.par <- matrix(nrow=length(unlist(start)), ncol=nind)

    y.list <- x.list <- opt <- list()
    yscores <- xscores  <- numeric(nr)
    f.min <- cor <- numeric(nind)

    for (i in 1:nind) {
        ## Do SVD on Y
        Yr = Y[ind[[i]],]
        nr=nrow(Yr)
        if (subject.smooth) {
            smoother <- function(ind, x, r){
                SVD <- svd(x[ind,], nu=r, nv=r)
                SVD$u %*% (SVD$d[1:r] * t(SVD$v))
            }
            if (global) ## smooth subject data separately, then diagonalise
                Yr <- bdiag(tapply(1:nrow(Yr), list(subject), 
                                   smoother, Yr, subject.smooth))
            else
                Yr <- smoother(1:nrow(Yr), Yr, subject.smooth)
        }
        ## Calculate RYr
        if (is.empty.model(partial)){
            y.list[[i]] = Yr
        }
        else { ## N.B equiv to doing for each subject separately, bit inefficient
            G <- model.matrix(mt, mf[ind[[i]], , drop = FALSE])
            y.list[[i]] = lm.fit(G, as.matrix(Yr))$residuals
        }   
        ## Reduce the dimensionality using svd
        rank.RYr=rankMatrix(y.list[[i]])
        SVD = svd(y.list[[i]], nu=rank.RYr, nv=rank.RYr )
        D<-diag(SVD$d[1:rank.RYr],ncol=rank.RYr)
        U1<-SVD$u
        V<-SVD$v
        RYr=U1%*%D
        S11=crossprod(RYr)
        S11.eig <- eigen(S11)
        S11.val=suppressWarnings(sqrt(S11.eig$values))
        S11.val[is.na(S11.val)]=0
        S11.sqrt <- solve(S11.eig$vectors %*% (S11.val * t(S11.eig$vectors)))
        dat <- mf[ind[[i]], , drop = FALSE]
        if (!is.null(treatment)) {
            ##can ignore reference level when finding nonlin par
            id <- unclass(factor(treatment[ind[[i]]]))

            F <- class.ind(id)
            if (!is.null(ref)) {
                F <- F[,-ref]
                nlev <- max(id)
                ord <- numeric(nlev)  
                ord[ref] <- nlev
                ord[-ref] <- seq(nlev - 1)
                id <- ord[id]
            }
            if (!separate) id <- 1
        }
        else F <- id <- 1
        nm <- names(start)
        qy <- qr(RYr)  
        dy <- qy$rank
        obj.f=function(Kvector) {
            ## evaluate RFr
            par <- split(Kvector, group)
            dat[nm] <- lapply(par, "[", id)
            val <- eval(formula, dat)
            val[is.na(val)] <- 0
            if (is.empty.model(partial)){
                RFr <- val*F
            }
            else{
                RFr <- lm.fit(G, val*F)$residuals
            }
            ## find cor as in MASS:::cancor
            qx <- qr(RFr)
            dx <- qx$rank
            z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , 
                     drop = FALSE], dx, dy)
            log(1-(z$d[1])^2)
        }
        ## Minimize ln(1-largesteigenvalue(Cor))
        opt[[i]] <- suppressWarnings(optim(unlist(start), obj.f, method = method,
                                           lower = lower, upper = upper, ...))
        nonlin.par[,i]=opt[[i]]$par
        par <- split(nonlin.par[,i], group)
        dat[nm] <- lapply(par, "[", id)
        val <- eval(formula, dat)
        val[is.na(val)] <- 0        
        RFr <- val*F
        if (!is.empty.model(partial))
            RFr <- lm.fit(G, RFr)$residuals
        ## find cor as in MASS:::cancor
        qx <- qr(RFr)
        dx <- qx$rank
        z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , 
                 drop = FALSE], dx, dy)
        ycoef[,i]= V %*% backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE],
                                   z$v)[,CCA.roots]
        B= backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)[,CCA.roots]
        RFrB=as.matrix(RFr)%*%B
        yscores[ind[[i]]] = as.matrix(y.list[[i]])%*%ycoef[,i]
        cor[i] <- lm.fit(RFrB, yscores[ind[[i]]])$coef
        f.min[i]=Re(opt[[i]]$value)
        ## Return RF, i.e. X partialling out G
        xscores[ind[[i]]] = RFrB
        x.list[[i]]=RFr
        ## Recalc fitted value based on starting values to get alignment right
        xcoef[,i] <- as.vector(B)
        par <- split(unlist(start), group)
        dat[nm] <- lapply(par, "[", id)
        val <- eval(formula, dat)
        val[is.na(val)] <- 0                
        RFr <- val*F
        if (!is.empty.model(partial))
            RFr <- lm.fit(G, RFr)$residuals
        if (cor(xscores[ind[[i]]], rowSums(as.matrix(RFr))) < 0) {
            xcoef[,i] <- -xcoef[,i]
            xscores[ind[[i]]] <- -xscores[ind[[i]]]
            ycoef[,i] <- -ycoef[,i]
            yscores[ind[[i]]] <- -yscores[ind[[i]]]
        }
    }

    if (!is.null(treatment)) 
        treatment_c <- paste("", levels(treatment)[setdiff(seq(nlev), ref)])
    else treatment_c <- ""
    if (separate) rownames(nonlin.par) = paste(group, treatment_c, sep = "")
    else rownames(nonlin.par) = group

    # need to name xcoef by columns of X (covariate columns partialled out)
    #if (!is.empty.model(partial)) 
    #    rownames(xcoef)=c(paste('formula', treatment_c), colnames(G))
    #else rownames(xcoef)=paste('formula', treatment_c)

    ## split signatures if necessary
    if (global == TRUE){
        ycoef <- matrix(ycoef, nrow=length(f))
    }
    rownames(ycoef) <- freq.nm

    if (!is.null(subject)) {
        if (!global){
            colnames(nonlin.par) <- colnames(xcoef) <-
                names(cum.pct.explained) <- names(cor) <- names(opt) <-
                names(y.list) <- names(x.list) <- 
                paste('subject',levels(subject))
        }
        colnames(ycoef) <- paste('subject',levels(subject))
    }

    out <- list(call = match.call(),
                ycoef = ycoef, # 'signatures' rows = subject, cols = freq
                xcoef = xcoef, # linear parameters on RHS
                yscores = yscores,
                xscores = xscores,
                subject = subject,
                treatment = treatment,
                time = mf$time,
                ref = ref, # reference level(could be NULL?)
                nonlinear.parameters = nonlin.par,
                cor = cor,
                y = y.list,
                x = x.list,
                global.smooth = global.smooth,
                subject.smooth = subject.smooth,
                pct.explained = cum.pct.explained,
                opt = opt)
    class(out) <- "gslcca"
    out
}

