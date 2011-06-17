ESLCCA <- function (Y, # matrix of power spectra
                    time, # time points
                    formula = "Double Exponential", # use built-in for now, implement formula later
                    subject = NULL, # subject indicator (allow NULL => no subjects/treat the same?)
                    treatment = NULL, # treatment factor (allow NULL => single curve?)
                    ref = 1, # reference level (allow NULL => no control?)
                    separate = TRUE, #fit separate nonlinear parameters for each treatment?
                    partial = ~1,
                    data = NULL,
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
                        "Double Exponential" = list(K1 = 8.5, K2 = 9),
                        "Critical Exponential" = list(K1 = 8.5))
        formula <- switch(formula,
                          "Double Exponential" = ~ (abs(K1-K2)>10e-6)*(exp(-time/exp(K1)) - exp(-time/exp(K2))) +
                              (abs(K1-K2)<=10e-6)*(time*exp(-time/exp(K1))) ,
                          "Critical Exponential" = ~ time*exp(-time/exp(K1)))
    }

    vars <- all.vars(formula)
    if (!"time" %in% vars) stop("'formula' must be a function of 'time'")

    ## anything without starting values assumed to be a variable
    vars <- setdiff(c(vars, all.vars(partial)), c("time", names(start)))
    dummy <- reformulate(c("0", vars)) #don't need intercept for model frame
    formula <- as.expression(formula[[length(formula)]])

    ## Collate data to ensure equal length and deal with NAs across all variables
    ## get model frame including Y, time, subject & treatment if specified, omit NAs
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("Y", "subject", "treatment", "time", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$formula <- dummy
    mf$na.action <- na.omit
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ## split up again for convenience
    Y <- mf$`(Y)`
    f <- suppressWarnings(as.numeric(colnames(Y)))
    if (any(is.na(f))) f <- seq_along(f)

    subject <- as.factor(mf$`(subject)`)
    if (is.null(subject)) subject <- 1 #is this enough?

    if (!is.null(mf$`(treatment)`)) {
        treatment <- as.factor(treatment)
        if (is.numeric(ref)) ref <- levels(treatment)[ref]
        include <- treatment != ref #can ignore reference level in computations
        id <- unclass(factor(treatment[include]))
        ntrt <- nlevels(treatment) - 1
    }
    else ntrt <- 1

    if (is.null(mf$`(time)`)) stop("'time' must be specified and have length equal to", nrow(Y), "\n")
    names(mf)[match("(time)", names(mf))] <- "time"

    ## leave only covariates in mf
    mf <- mf[, !(names(mf) %in% c("(Y)", "(time)", "(subject)", "(treatment)")), drop = FALSE]
    mt <- terms(partial)
    if (!is.empty.model(partial)) {
        G <- model.matrix(mt, mf[1,]) #temp value to find number of parameters associated with covariates
        np <- ncol(G)
    }
    else np <- 0
    nind=nlevels(subject)
    xcoef=matrix(, nr = ntrt + np, ncol = nind)

    ## Pre-smoothing of the complete data set for artefacts removal
    if (global.smooth) {
        if (global.smooth == TRUE) {
            Ytemp= scale(Y,scale=F)
            eig=cbind(eigen(crossprod(Ytemp))$values)
            nroots=1:length(eig)
            global.smooth=max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]) + 1
        }
        SVD = svd( Y, nu=global.smooth, nv=global.smooth )
        Y = SVD$u%*%diag(SVD$d[1:global.smooth],ncol=global.smooth)%*%t(SVD$v)
    }

    ycoef=matrix(nrow=nind,ncol=length(f))
    n <- nrow(Y)
    yscores = numeric(n)

    if (subject.smooth == TRUE) {
        ## Automatically select number of roots
        subject.smooth <- numeric(nind)
        for (i in seq_len(nind)) {
            Ytemp = scale(Y[subject==levels(subject)[i],],scale=F)
            eig = eigen(crossprod(Ytemp), only.values = TRUE)
            nroots = seq_along(eig)
            subject.smooth[i] = max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98])
        }
        subject.smooth = max(subject.smooth)+1
    }
    else if (!is.numeric(subject.smooth)) {
        ## No smoothing case, so r will be equal to the rank of the data matrix
        subject.smooth = min(c(min(table(subject)),length(f)))
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

    CCA.roots = 1

    nonlin.par=matrix(nrow=nind,ncol=length(unlist(start)))

    f.min=numeric(nind)
    y.list=list()
    x.list=list()
    opt <- list()

    corr=numeric(nind)
    xscores=numeric(n)

    for (i in 1:nind) {
        ## Do SVD on Y
        Indsubject = subject==levels(subject)[i]
        Yr = Y[Indsubject,]
        nr=nrow(Yr)
        SVD = svd( Yr, nu=subject.smooth, nv=subject.smooth )
        Yr = SVD$u%*%tcrossprod(x=diag(SVD$d[1:subject.smooth]),y=SVD$v)
        ## Calculate the idempotent matrix R
        if (is.empty.model(partial)){
            R = diag(nr)
        }
        else{
            G <- model.matrix(mt, subset(mf, Indsubject))
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
        S11.sqrt <- solve(S11.eig$vectors %*% tcrossprod(diag(S11.val,ncol=length(S11.val)),S11.eig$vectors) )
        if (!is.null(treatment)) {
            include <-  treatment[Indsubject] != ref #can ignore reference level when finding nonlin par
            id <- unclass(factor(treatment[Indsubject][include]))
            R <- R[,include]
            dat <- subset(mf, Indsubject & include)
            F <- class.ind(id)
            if (!separate) id <- 1
        }
        else {
            dat <- subset(mf, Indsubject)
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
        ## Minimize ln(1-largesteigenvalue(Corr))
        opt[[i]] <- optim(unlist(start), obj.f, ...)
        nonlin.par[i,]=opt[[i]]$par
        par <- split(nonlin.par[i,], group)
        dat[nm] <- lapply(par, "[", id)
        val <- eval(formula, dat)
        Fr <- val*F
        RFr <- R %*% Fr
        S21 <- t(RFr) %*% RYr
        S22.inv.S21 <- solve(crossprod(RFr)) %*% S21
        Eigvectors=eigen(S11.sqrt %*% t(S21) %*% S22.inv.S21 %*% S11.sqrt, symmetric = TRUE)$vectors
        ycoef[i,]= Re(V%*%S11.sqrt%*%Eigvectors[,CCA.roots])
        B= Re(S22.inv.S21%*%S11.sqrt%*%Eigvectors[,CCA.roots])
        RFrB=RFr%*%B
        Norm <- sqrt(colSums(RFrB^2))
        B <- B/matrix(Norm,nrow(B),ncol=CCA.roots,byrow=T)
        RFrB=RFr%*%B
        y.list[[i]]=Yr
        yscores[Indsubject] = y.list[[i]]%*%ycoef[i,]
        f.min[i]=Re(opt[[i]]$value)
        ## Save intercept and covariates linear parameters.
        if (length(include) > 1) x.list[[i]]=rbind(matrix(0,nrow=sum(!include),ncol=ncol(Fr)),Fr)
        else x.list[[i]]=Fr
        xscores[Indsubject] = x.list[[i]]%*%B
        corr[i] <- cor(yscores[Indsubject], xscores[Indsubject])
        if (!is.empty.model(partial)) {
            x.list[[i]]=cbind(x.list[[i]],G)
            lin <- lm.fit(G, yscores[Indsubject]/corr[i] - xscores[Indsubject])
            xcoef[,i] <- c(B, lin$coef)
            xscores[Indsubject] = xscores[Indsubject] + fitted(lin)
        }
        else xcoef[,i] <- B
        if (cor(yscores[Indsubject], rowSums(y.list[[i]])) < 0) {
            xcoef[,i] <- -xcoef[,i]
            xscores[Indsubject] <- -xscores[Indsubject]
            ycoef[i,] <- -ycoef[i,]
            yscores[Indsubject] <- -yscores[Indsubject]
        }
    }

    if (!is.null(treatment)) treatment_c <- paste("", setdiff(levels(treatment), ref))
    else treatment_c <- ""
    rownames(nonlin.par)= paste('subject',levels(subject))
    colnames(nonlin.par)= paste(rep(paste('Parameter',1:nPar),rep(reps,nPar)),rep(ifelse(separate, treatment_c, ""),nPar), sep = "")

    if (!is.null(subject)) colnames(xcoef)=as.character(paste('subject',levels(subject)))
    if (!is.empty.model(partial)) rownames(xcoef)=c(paste('formula', treatment_c), colnames(G))
    else rownames(xcoef)=paste('formula', treatment_c)

    ycoef <- as.data.frame(ycoef)
    colnames(ycoef) <- f
    if (!is.null(subject)) rownames(ycoef) <- levels(subject)

    R.square = matrix(1-exp(f.min),ncol=1, dimnames=(list(rownames(nonlin.par),'R^2')))

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
                corr = corr,
                y = y.list,
                x = x.list,
                global.roots = global.smooth,
                subject.roots = subject.smooth,
                opt = opt)
    class(out) <- "ESLCCA"
    out
}

