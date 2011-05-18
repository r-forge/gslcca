ESLCCA <- function (Y, # matrix of power spectra
                    time, # time points
                    formula = 'DoubleExponential', # use built-in for now, implement formula later
                    subject = NULL, # subject indicator (allow NULL => no subjects/treat the same?)
                    treatment = NULL, # treatment factor (allow NULL => single curve?)
                    ref = 1, # reference level (allow NULL => no control?)
                    na.action,
                    data.roots=0,subject.roots=-1,con=TRUE,
                    separate=TRUE,
                    start = NULL,
                    ...) { # allow arguments to be passed to optim

# Upload required R libraries
library(MASS)

# Prepare the data
if (missing(na.action)) na.action <- getOption("na.action")
Y <- match.fun(na.action)(Y)
f <- suppressWarnings(as.numeric(colnames(Y)))
if (any(is.na(f))) f <- seq_along(f)

n <- nrow(Y)
if (!(length(subject) %in% c(0, n))) stop("length of 'subject' must equal nrow(Y)")
if (!(length(treatment) %in% c(0, n))) stop("length of 'treatment' must equal nrow(Y)")
if (missing(time) || length(time) != n) stop("'time' must be specified and have length equal to nrow(Y)")

if (length(attr(Y, "na.action"))){
    if (!missing(subject)) subject <- subject[-attr(Y, "na.action")]
    if (!missing(treatment)) treatment <- treatment[-attr(Y, "na.action")]
    time <- time[-attr(Y, "na.action")]
}

##### assume treatment non-NULL for now #######
treatment <- as.factor(treatment)
if (is.numeric(ref)) ref <- levels(treatment)[ref]
include <- treatment != ref #can ignore reference level in computations
id <- unclass(factor(treatment[include]))

subject <- as.factor(subject)

# Pre-smoothing of the complete data set for artifacts removal
if (data.roots!=0) {
    if  (data.roots==-1) {
         Ytemp= scale(Y,scale=F)
         eig=cbind(eigen(crossprod(Ytemp))$values)
         nroots=1:length(eig)
         data.roots=max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]) + 1
    }
    SVD = svd( Y, nu=data.roots, nv=data.roots )
    Y = SVD$u%*%diag(SVD$d[1:data.roots],ncol=data.roots)%*%t(SVD$v)
}

Y=Y-apply(Y,1,min)
nind=nlevels(subject)
ntreat=nlevels(treatment)
Gamma=matrix(nrow=nind,ncol=length(f))
observed.val = c()
if (separate==TRUE) {
    nvar=ntreat-1
    size = function() {n_vec}
} else {
    nvar=1
    size = function() {nr-nv}
}

if (subject.roots<=0) {
    if (subject.roots==0) {
        # No smoothing case, so r will be equal to the rank of the data matrix
        subject.roots = min(c(min(table(subject)),length(f)))
    } else {
        # Default smoothing
        subject.roots = c()
        eig = c()
        for (i in 1:nind) {
             Ytemp = scale(Y[subject==levels(subject)[i],],scale=F)
             eig = cbind(eigen(crossprod(Ytemp))$values)
             nroots = 1:length(eig)
             subject.roots = c(subject.roots,max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]))
        }
        subject.roots = max(subject.roots)+1
    }
}

# In the case that the curve is not specified by the user, then the default is
# the double exponential (with positive rate parameters) and the default initial
# values are following
if (is.character(formula)) {
    start <- switch(formula,
                    "DoubleExponential" = list(K1 = rep(8.5, nvar), K2 = rep(9, nvar)),
                    "CriticalExponential" = list(K1 = rep(8.5, nvar)))
}

if (is.null(start)) {
    stop('No initial values have been specified for the nonlinear parameters.')
}
else {
    start <- mapply(function(start, label, nvar) {
        if (length(start) == 1) rep(start, nvar)
        else if (length(start) == nvar) start
        else stop("there should be 1 or", nvar, "starting values for parameter", label, "\n")
    }, start = start, label = names(start), MoreArgs = list(nvar = nvar), SIMPLIFY = FALSE)
}
nPar <- length(start)
group <- gl(nPar, nvar, labels = names(start))

if (is.character(formula)) {
    formula <- switch(formula,
                      "DoubleExponential" = expression( (abs(K1-K2)>10e-6)*(exp(-time/exp(K1)) - exp(-time/exp(K2))) +
                          (abs(K1-K2)<=10e-6)*(time*exp(-time/exp(K1))) ),
                      "CriticalExponential" = expression(time*exp(-time/exp(K1))))
}
else formula <- as.expression(formula[[length(formula)]])
CCA.roots = 1

nonlin.par=matrix(nrow=nind,ncol=length(unlist(start)))

Theta.m=c()
linear.par = C = G = c()
f.min=c()
y.list=list()
x.list=list()
opt <- list()

corr_mat=matrix(ncol=nind,nrow=nlevels(treatment))
fitted.val=c()

if (!is.numeric(con)) {
      if (con==TRUE) {
          case = 1
      } else {
          case = 2
      }
  } else {
      case = 3
}

for (i in 1:nind) {
     # Do PCA on Y
     Indsubject = subject==levels(subject)[i]
     Yr = Y[Indsubject,]
     nr=nrow(Yr)
     SVD = svd( Yr, nu=subject.roots, nv=subject.roots )
     Yr = SVD$u%*%tcrossprod(x=diag(SVD$d[1:subject.roots]),y=SVD$v)
     # Calculate the idempotent matrix R
     I = rep(1,nr)
     if (!is.numeric(con)) {
         if (case==1) {
             G=I
             R = diag(nr) - matrix(1/nr, nr, nr)
         } else {
             R = diag(I)
         }
     } else {
         G=con[Indsubject,]
         R = diag(I)-G%*%solve(crossprod(G))%*%t(G)
     }
     RYr=R%*%Yr
     # Reduce the dimensionality using svd
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
     include <-  treatment[Indsubject] != ref #can ignore reference level when finding nonlin par
     id <- unclass(factor(treatment[Indsubject][include]))
     R2=R
     R <- R[,include]
     dat <- data.frame(time = time[Indsubject][include])
     nm <- names(start)
     env <- environment(dat)
     F <- class.ind(id)
     obj.f=function(Kvector) {
         par <- split(Kvector, group)
         dat[nm] <- lapply(par, "[", id)
         val <- eval(formula, dat)
         RFr <- R %*% (val*F)
         C <- chol(crossprod(RFr)) # cholesky decomposition of RFr %*% t(RFr)
         z <- backsolve(C, crossprod(RFr, RYr) %*% S11.sqrt, transpose = TRUE) # solve t(C) %*% S21 %*% S11.sqrt for z
         log(1-eigen(crossprod(z), symmetric = TRUE, only.values = TRUE)$values[1])
     }
     # Minimize ln(1-largesteigenvalue(Corr))
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
     optEigvalue=Eigvectors[,CCA.roots]
     Gamma[i,]= Re(V%*%S11.sqrt%*%optEigvalue)
     Theta= Re(S22.inv.S21%*%Eigvectors)
     RFrTheta=RFr%*%Theta
     Norm <- colSums(RFrTheta^2)
     Theta=Theta/sqrt(matrix(Norm,nrow(Theta),ncol=ncol(Theta),byrow=T))
     RFrTheta=RFr%*%Theta[,CCA.roots]
     hRFrTheta <- RFrTheta[include][-(1:10)]
     cond=hRFrTheta[abs(hRFrTheta)==max(abs(hRFrTheta))]
    # if   (cond[1]<0) {Gamma[i,]=-Gamma[i,]; RFrTheta=-RFrTheta}
     y.list[[i]]=tcrossprod(RYr,V)
     observed.val=c( observed.val, y.list[[i]]%*%Gamma[i,] )
     fitted.val=c(fitted.val,RFrTheta)
     Theta.m=cbind(Theta.m,Theta[,1])
     f.min=c(f.min,Re(opt[[i]]$value))
     x.list[[i]]=cbind(RFr,G)
     # Save intercept and covariates linear parameters.
     Y.tilde=ginv(R2)%*%Yr%*%as.matrix(Gamma[i,])
     F.tilde=c(rep(0,sum(include==0)),Fr%*%Theta[,1])
     if (case != 2) { C=solve(crossprod(G))%*%t(G)%*%(Y.tilde-F.tilde) }
     linear.par=cbind(linear.par,C)
}

l_ind=levels(subject)
mydata.pca <- princomp(t(Gamma))
loading <- mydata.pca$loadings[,1]
sign_sign= loading>0
sl_ind=l_ind[sign_sign]
if (sum(sign_sign)<=sum(sign_sign==F)) {
     Gamma[sign_sign,]=-Gamma[sign_sign,]
     fitted.val[subject %in% sl_ind] = - fitted.val[subject %in% sl_ind]
     observed.val[subject %in% sl_ind] = - observed.val[subject %in% sl_ind]
}  else {
sl_ind= l_ind[sign_sign==F]
     Gamma[sign_sign==F,]=-Gamma[sign_sign==F,]
     fitted.val[subject %in% sl_ind] = - fitted.val[subject %in% sl_ind]
     observed.val[subject %in% sl_ind] = - observed.val[subject %in% sl_ind]
}

treatment_c <- setdiff(levels(treatment), ref)
rownames(nonlin.par)=as.character(paste('subject',levels(subject)))
if (separate==T) {
colnames(nonlin.par)= paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)),rep(treatment_c,nPar))
} else {
colnames(nonlin.par)=paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)))
}

if (case != 1) {
    rownames(linear.par) = colnames(con)
} else {
    rownames(linear.par) = 'Intercept'
}
colnames(Theta.m)=as.character(paste('subject',levels(subject)))
rownames(Theta.m)=treatment_c
xcoef = rbind(Theta.m, linear.par)

ycoef <- as.data.frame(Gamma)
colnames(ycoef) <- f
rownames(ycoef) <- levels(subject)

R.square = matrix(1-exp(f.min),ncol=1, dimnames=(list(rownames(nonlin.par),'R^2')))

out <- list(call = match.call(),
            ycoef = ycoef, # 'signatures' rownames = subject, colnames = freq
            xcoef = xcoef, # yet to implement, should be linear parameters on RHS; rownames - trt levels + cov names + Intercept
            scores = data.frame(subject = subject,
                                treatment = treatment,# could be NULL?
                                time = time,
                                xscores = fitted.val,
                                yscores = observed.val),
            ref = ref, # could be NULL?
            nonlinear.parameters = as.data.frame(nonlin.par),
            R.square = R.square,
            y = y.list,
            x = x.list,
            global.roots = data.roots,
            subject.roots = subject.roots,
            opt = opt)
class(out) <- "ESLCCA"
out
}

