ESLCCA <- function (Y, # matrix of power spectra
                    formula = "Double Exponential", # use built-in for now, implement formula later
                    time, # time points
                    subject = NULL, # subject indicator (allow NULL => no subjects/treat the same?)
                    treatment = NULL, # treatment factor (allow NULL => single curve?)
                    ref = 1, # reference level (allow NULL => no control?)
                    separate = TRUE, #fit separate nonlinear parameters for each treatment?
                    partial = ~1,
                    data = NULL,
                    subset = NULL,
                    global.smooth = FALSE,
                    subject.smooth = TRUE,
                    start = NULL,
                    ...) { # allow arguments to be passed to optim
    ## Formula must be a function of "time", possibly other var and parameters
    if (is.character(formula)){
        models <- c("Double Exponential", "Critical Exponential")
        formula <- models[agrep(formula, models)]
        if (!length(formula)) stop("\"formula\" not recognised, only \"Double Exponential\" or \"Critical Exponential\" \n",
                                   "can be specified as character strings")
        start <- switch(formula, #replicate later for multiple treatments
                        "Double Exponential" = list(K1 = 9, K2 = 8.5), #K1 > K2 so increases from ref
                        "Critical Exponential" = list(K1 = 8.5))
        formula <- switch(formula,
                          "Double Exponential" = ~ (abs(K2-K1)>10e-6)*(exp(-time/exp(K1)) - exp(-time/exp(K2))) +
                              (abs(K2-K1)<=10e-6)*(time*exp(-time/exp(K1))) ,
                          "Critical Exponential" = ~ time*exp(-time/exp(K1)))
    }

    vars <- all.vars(formula)
    if (!"time" %in% vars) stop("'formula' must be a function of 'time'")

    ## anything without starting values assumed to be a variable
    vars <- setdiff(c(vars, all.vars(partial)), c("time", names(start)))
    dummy <- reformulate(c("0", vars), response = "Y") #don't need intercept for model frame
    formula <- as.expression(formula[[length(formula)]])

    ## Collate data to ensure equal length and deal with NAs across all variables
    ## get model frame including Y, time, subject & treatment if specified, omit NAs
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("subject", "treatment", "time", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- dummy
    mf$na.action <- na.omit
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## split up again for convenience
    Y <- mf$Y
    nr <- nrow(Y)
    f <- suppressWarnings(as.numeric(colnames(Y)))
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
        if (is.numeric(ref)) ref <- levels(treatment)[ref]
        include <- treatment != ref #can ignore reference level in computations
        id <- unclass(factor(treatment[include]))
        ntrt <- nlevels(treatment) - 1
    }
    else ntrt <- 1

    if (is.null(mf$`(time)`)) stop("'time' must be specified and have length equal to", nr, "\n")
    names(mf)[match("(time)", names(mf))] <- "time"

    ## take out extras from mf
    mf <- mf[, !(names(mf) %in% c("Y", "(subject)", "(treatment)")), drop = FALSE]

    ## Pre-smoothing of the complete data set for artefacts removal
    if (global.smooth) {
        if (global.smooth == TRUE) {
            Ytemp= scale(Y,scale=F)
            eig=eigen(crossprod(Ytemp), only.values = TRUE)$values
            nroots=seq_along(eig)
            global.smooth=max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]) + 1
        }
        SVD = svd( Y, nu=global.smooth, nv=global.smooth )
        Y = SVD$u%*% (SVD$d[1:global.smooth] * t(SVD$v))
    }

    if (isTRUE(subject.smooth)) {
        ## Automatically select number of roots
        subject.smooth <- numeric(nind)
        for (i in seq_len(nind)) {
            Ytemp = scale(Y[ind[[i]],],scale=F)
            eig = eigen(crossprod(Ytemp), only.values = TRUE)$values
            nroots = seq_along(eig)
            subject.smooth[i] = max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98])
        }
        subject.smooth = max(subject.smooth)+1
    }

    reps <- ifelse(separate, ntrt, 1)
    if (is.null(start)) {
        stop('No initial values have been specified for the nonlinear parameters.')
    }
    else {
        start <- mapply(function(start, label, reps) {
            if (length(start) == 1) rep(start, reps)
            else if (length(start) == reps) start
            else stop("there should be 1 or", reps, "starting values for parameter", label, "\n")
        }, start = start, label = names(start), MoreArgs = list(reps = reps), SIMPLIFY = FALSE)
    }
    nPar <- length(start)
    group <- gl(nPar, reps, labels = names(start))

    ## evaluate RHS for first row of data
    par <- as.list(rep(NA, length(start)))
    names(par) <- names(start)
    dat <- do.call("cbind", list(mf[1,, drop = FALSE], par))
    val <- eval(formula, dat)
    mt <- terms(partial)
    if (!is.empty.model(partial)) {
        G <- model.matrix(mt, mf[1,]) #temp value to find number of parameters associated with covariates
        np <- ncol(G)
    }
    else np <- 0
    CCA.roots = 1

    nonlin.par <- matrix(nrow=length(unlist(start)), ncol=nind)
    ycoef <- matrix(nrow=length(f), ncol=nind)
    xcoef <- matrix(nrow=length(val)*ntrt + np, ncol=nind)

    y.list <- x.list <- opt <- list()
    yscores <- xscores  <- numeric(nr)
    f.min <- cor <- numeric(nind)

    for (i in 1:nind) {
        ## Do SVD on Y
        Yr = Y[ind[[i]],]
        nr=nrow(Yr)
        if (subject.smooth) {
            SVD = svd( Yr, nu=subject.smooth, nv=subject.smooth )
            Yr = SVD$u %*% (SVD$d[1:subject.smooth] * t(SVD$v))
        }
        ## Calculate the idempotent matrix R
        if (is.empty.model(partial)){
            R = diag(nr)
        }
        else{
            G <- model.matrix(mt, mf[ind[[i]], , drop = FALSE])
            R = diag(nr)-G%*%solve(crossprod(G))%*%t(G)
        }
        RYr=R%*%Yr ## better to compute directly as resids?
        ## Reduce the dimensionality using svd
        rank.RYr=qr(RYr)$rank
        SVD = svd( RYr, nu=rank.RYr, nv=rank.RYr )
        D<-diag(SVD$d[1:rank.RYr],ncol=rank.RYr)
        U1<-SVD$u
        V<-SVD$v
        RYr=U1%*%D
        S11=crossprod(RYr)
        S11.eig <- eigen(S11)
        S11.val=suppressWarnings(sqrt(S11.eig$values))
        S11.val[is.na(S11.val)]=0
        S11.sqrt <- solve(S11.eig$vectors %*% (S11.val * t(S11.eig$vectors)))
        if (!is.null(treatment)) {
            include <-  treatment[ind[[i]]] != ref #can ignore reference level when finding nonlin par
            id <- unclass(factor(treatment[ind[[i]]][include]))
            R <- R[,include]
            dat <- mf[ind[[i]][include], , drop = FALSE]
            F <- class.ind(id)
            if (!separate) id <- 1
        }
        else {
            dat <- mf[ind[[i]], , drop = FALSE]
            F <- id <- 1
            include <- TRUE
        }
        nm <- names(start)
        obj.f=function(Kvector) {
            par <- split(Kvector, group)
            dat[nm] <- lapply(par, "[", id)
            val <- eval(formula, dat)
            RFr <- R %*% (val*F)
            C <- chol(crossprod(RFr)) # cholesky decomposition of RFr %*% t(RFr)
            z <- backsolve(C, crossprod(RFr, RYr) %*% S11.sqrt, transpose = TRUE) # solve t(C) %*% S21 %*% S11.sqrt for z
            log(1-eigen(crossprod(z), symmetric = TRUE, only.values = TRUE)$values[1])
        }
        ## Minimize ln(1-largesteigenvalue(Cor))
        opt[[i]] <- optim(unlist(start), obj.f, ...)
        nonlin.par[,i]=opt[[i]]$par
        par <- split(nonlin.par[,i], group)
        dat[nm] <- lapply(par, "[", id)
        val <- eval(formula, dat)
        Fr <- val*F
        RFr <- R %*% Fr
        S21 <- t(RFr) %*% RYr
        S22.inv.S21 <- solve(crossprod(RFr)) %*% S21
        Eigvectors=eigen(S11.sqrt %*% t(S21) %*% S22.inv.S21 %*% S11.sqrt, symmetric = TRUE)$vectors
        ycoef[,i]= Re(V%*%S11.sqrt%*%Eigvectors[,CCA.roots])
        B= Re(S22.inv.S21%*%S11.sqrt%*%Eigvectors[,CCA.roots])
        RFrB=RFr%*%B
        Norm <- sqrt(colSums(RFrB^2))
        B <- B/matrix(Norm,nrow(B),ncol=CCA.roots,byrow=T)
        RFrB=RFr%*%B
        y.list[[i]]=Yr
        yscores[ind[[i]]] = y.list[[i]]%*%ycoef[,i]
        cor[i] <- lm(yscores[ind[[i]]] ~ 0 + RFrB)$coef
        f.min[i]=Re(opt[[i]]$value)
        ## Save intercept and covariates linear parameters.
        if (length(include) > 1) x.list[[i]]=rbind(matrix(0,nrow=sum(!include),ncol=ncol(Fr)),Fr)
        else x.list[[i]]=Fr
        xscores[ind[[i]]] = x.list[[i]]%*%B
        if (!is.empty.model(partial)) {
            x.list[[i]]=cbind(x.list[[i]],G)
            lin <- lm.fit(G, yscores[ind[[i]]] - xscores[ind[[i]]])
            xcoef[,i] <- c(B, lin$coef)
            xscores[ind[[i]]] = xscores[ind[[i]]] + fitted(lin)
        }
        else xcoef[,i] <- B
        par <- split(unlist(start), group)
        dat[nm] <- lapply(par, "[", id)
        val <- eval(formula, dat)
        Fr <- val*F
        RFr <- R %*% Fr
        if (cor(xscores[ind[[i]]], rowSums(RFr)) < 0) {
            xcoef[,i] <- -xcoef[,i]
            xscores[ind[[i]]] <- -xscores[ind[[i]]]
            ycoef[,i] <- -ycoef[,i]
            yscores[ind[[i]]] <- -yscores[ind[[i]]]
        }
    }

    if (!is.null(treatment)) treatment_c <- paste("", setdiff(levels(treatment), ref))
    else treatment_c <- ""
    rownames(nonlin.par)= paste(rep(names(start),rep(reps,nPar)),rep(ifelse(separate, treatment_c, ""),nPar), sep = "")

    if (!is.empty.model(partial)) rownames(xcoef)=c(paste('formula', treatment_c), colnames(G))
    else rownames(xcoef)=paste('formula', treatment_c)

    ycoef <- as.data.frame(ycoef)
    rownames(ycoef) <- f

    if (!is.null(subject))
        colnames(nonlin.par) <- colnames(xcoef) <- colnames(ycoef) <- paste('subject',levels(subject))

    R.square = matrix(1-exp(f.min),ncol=1, dimnames=(list(colnames(nonlin.par),'R^2')))

    out <- list(call = match.call(),
                ycoef = ycoef, # 'signatures' rownames = subject, colnames = freq
                xcoef = xcoef, # yet to implement, should be linear parameters on RHS; rownames - trt levels + cov names + Intercept
                yscores = yscores,
                xscores = xscores,
                subject = subject,
                treatment = treatment,
                time = mf$time,
                ref = ref, # reference level, not implemented yet (could be NULL?)
                nonlinear.parameters = as.data.frame(nonlin.par),
                R.square = R.square,
                cor = cor,
                y = y.list,
                x = x.list,
                global.roots = global.smooth,
                subject.roots = subject.smooth,
                opt = opt)
    class(out) <- "ESLCCA"
    out
}

