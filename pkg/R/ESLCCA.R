ESLCCA <- function (Y, # matrix of power spectra
                    subject = NULL, # subject indicator (allow NULL => no subjects/treat the same?)
                    treatment = NULL, # treatment factor (allow NULL => single curve?)
                    time, # time points
                    ref = 1, # reference level (allow NULL => no control?)
                    na.action,
                    Curve = 'DoubleExponential', # use built-in for now, implement formula later
                    PreSmoothingRoots=0,r=-1,con=0,separate=T,path='c:\\',comp_name,
                    start = NULL) {

# Upload required R libraries
library(MASS)

# Prepare the data
if (missing(na.action)) na.action <- getOption(na.action)
Y <- na.action(Y)
f <- as.numeric(colnames(Y))
if (any(is.na(f))) f <- seq_along(f)

n <- nrow(Y)
if (!(length(subject) %in% c(0, n))) stop("length of 'subject' must equal nrow(Y)")
if (!(length(treatment) %in% c(0, n))) stop("length of 'treatment' must equal nrow(Y)")
if (length(time) != n) stop("'time' must be specified and have length equal to nrow(Y)")

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

# Pre-smoothing of the complete data set for artifacts removal
if (PreSmoothingRoots!=0) {
    if  (PreSmoothingRoots==-1) {
         Ytemp= scale(Y,scale=F)
         eig=cbind(eigen(t(Ytemp)%*%Ytemp)$values)
         nroots=1:length(eig)
         PreSmoothingRoots=max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]) + 1
    }
    SVD = svd( Y, nu=PreSmoothingRoots, nv=PreSmoothingRoots )
    Y = SVD$u%*%diag(SVD$d[1:PreSmoothingRoots])%*%t(SVD$v)
}

Y=Y-apply(Y,1,min)
nind=nlevels(subject)
ntreat=nlevels(treatment)
Gamma=matrix(nrow=nind,ncol=length(f))
Fitted_values=treatment_ind=Time_ind= c()
if (separate==T) {
    nvar=ntreat-1
    size = function() {n_vec}
} else {
    nvar=1
    size = function() {nr-nv}
}

if (r<=0) {
    if (r==0) {
        # No smoothing case, so r will be equal to the rank of the data matrix
        r=min(c(min(table(subject)),length(f)))
    } else {
        # Default smoothing
        r=c()
        eig=c()
        for (i in 1:nind) {
             Ytemp= scale(Y[subject==levels(subject)[i],],scale=F)
             eig=cbind(eigen(t(Ytemp)%*%Ytemp)$values)
             nroots=1:length(eig)
             r=c(r,max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]))
        }
        r=max(r)+1
    }
}

# In the case that the curve is not specified by the user, then the default is
# the double exponential (with positive rate parameters) and the default initial
# values are following
if (is.null(start)) {
    if (!is.function(Curve)) {
        if (Curve == 'DoubleExponential') start <- list(K1 = rep(8.5, nvar), K2 = rep(9, nvar))
        if (Curve == 'CriticalExponential') start <- list(K1 = rep(8.5, nvar))
    }else{
        stop('No initial values have been specified for the nonlinear parameters.')
    }
}
else {
    start <- mapply(function(start, label, nvar) {
        if (length(start) == 1) rep(start, nvar)
        else if (length(start) == nvar) start
        else stop("there should be 1 or", nvar, "starting values for parameter", label, "\n")
    }, start = start, label = names(start), MoreArgs = list(nvar = nvar))
}
group <- gl(length(start), nvar)

if (!is.function(Curve)) {
    if (Curve == 'DoubleExponential') {
        Curve <- expression(exp(-time/exp(K1))-exp(-time/exp(K2)))
    }
    else if (Curve == 'CriticalExponential') {
        Curve <- expression(time*exp(-time/exp(K1)))
    }
}

root=1

nlPAR=matrix(nrow=nind,ncol=length(initial.K))

Theta_1=c()
f_min=c()
RYr_list=list()
RFr_list=list()

subject_ind=c()
corr_mat=matrix(ncol=nind,nrow=nlevels(treatment))
observed_val=c()
for (i in 1:nind) {
     t=c()
     n_vec=c()
     # Do PCA on Y
     Indsubject = subject==levels(subject)[i]
     subject_ind = c(subject_ind,as.character(subject)[Indsubject])
     IndV = Indsubject==T & treatment=='Control'
     nv = sum(IndV); tv = Time[IndV]
     Yr = Y[IndV,]
     treatment_c=as.character(levels(treatment))
     treatment_c=treatment_c[treatment_c!='Control']
     char.treatment=as.character(treatment)
     treatment_ind= c(treatment_ind,char.treatment[IndV])
     F0=rep(rep(0,length(treatment_c)),nv)
     for (j in 1:length(treatment_c))  {
     IndL = Indsubject==T & treatment==treatment_c[j]
     nl = sum(IndL); tl = Time[IndL]
     Yr = rbind(Yr,Y[IndL,])
     n_vec=c(n_vec,nl)
     t=c(t,tl)
     treatment_ind= c(treatment_ind,char.treatment[IndL])
     new_F=rep(0,length(treatment_c))
     new_F[j]=1
     F0=c(F0,rep(new_F,nl))
     }
     Time_ind=c(Time_ind,c(tv,t))
     nr=nrow(Yr)
     SVD = svd( Yr, nu=r, nv=r )
     Yr = SVD$u%*%diag(SVD$d[1:r])%*%t(SVD$v)
     # Calculate the idempotent matrix R
     I = rep(1,nr)
     if (con==0) {
         R = diag(I)-I%*%solve(t(I)%*%I)%*%t(I)
     } else {
         R = diag(I)
     }
     RYr=R%*%Yr
     # Reduce the dimensionality using svd
     rank_RYr=qr(RYr)$rank
     SVD = svd( RYr, nu=rank_RYr, nv=rank_RYr )
     D<-diag(SVD$d[1:rank_RYr])
     U1<-SVD$u
     V<-SVD$v
     RYr=U1%*%D
     # Define F
     F0=matrix(F0,nrow=nr,byrow=T)
     Fnl=function(Parameters) {
         Fb=F0
         Fb[F0==1]=Curve(t,Parameters)
         Fb
         }
     S11=t(RYr)%*%RYr
     S11.eig <- eigen(S11)
     S11.val=suppressWarnings(sqrt(S11.eig$values))
     S11.val[is.na(S11.val)]=0
     S11.sqrt <- solve(S11.eig$vectors %*% diag(S11.val) %*% t(S11.eig$vectors))
     nPar=length(initial.K)/length(size())  # computes the number of parameters in 'Curve'
     temp1=rep(size(),nPar)
     R <- R[, include]
     obj.f=function(Kvector) {
         par <- split(Kvector, group)
         dat <- lapply(par, "[", id)
         val <- eval(Curve, dat) #finding time?
         RFr <- rowsum(val * R, treatment)
         RFr = R%*%Fnl(matrix(rep(Kvector,temp1),ncol=nPar))
         S12=t(RYr)%*%RFr
         log(1-eigen(S11.sqrt%*%S12%*%ginv(t(RFr)%*%RFr)%*%t(S12)%*%S11.sqrt)$values[1])
     }
     # Minimize ln(1-largesteigenvalue(Corr))
     nlPAR[i,]=optim(initial.K,obj.f,control=list(warnOnly = T))$par
     RFr = R%*%Fnl(matrix(rep(nlPAR[i,],temp1),ncol=nPar))
     S12=t(RYr)%*%RFr
     S22.inv.S21=ginv(t(RFr)%*%RFr)%*%t(S12)
     Eigvectors=eigen(S11.sqrt%*%S12%*%S22.inv.S21%*%S11.sqrt)$vectors
     optEigvalue=Eigvectors[,root]
     Gamma[i,]= Re(V%*%S11.sqrt%*%optEigvalue)
     Gamma= Gamma
     Theta= Re(S22.inv.S21%*%Eigvectors)
     RFrTheta=RFr%*%Theta
     Norm=apply(RFrTheta^2,2,sum)
     Theta=Theta/sqrt(matrix(Norm,nrow(Theta),ncol=ncol(Theta),byrow=T))
     RFrTheta=RFr%*%Theta[,root]
     hRFrTheta=RFrTheta[-c(1:(length(t)+10))]
     cond=hRFrTheta[abs(hRFrTheta)==max(abs(hRFrTheta))]
     if   (cond[1]<0) {Gamma[i,]=-Gamma[i,]; RFrTheta=-RFrTheta}
     Fit_values=RYr%*%t(V)%*%Gamma[i,]
     Fitted_values=c(Fitted_values,Fit_values)
     observed_val=c(observed_val,RFrTheta)
     Theta_1=cbind(Theta_1,Theta[,1])
     f_min=c(f_min,Re(obj.f(nlPAR[i,])))
     RYr_list[[i]]=RYr%*%t(V)
     RFr_list[[i]]=RFr
}

l_ind=levels(subject)
mydata.pca <- princomp(t(Gamma))
loading <- mydata.pca$loadings[,1]
sign_sign= loading>0
sl_ind=l_ind[sign_sign]
if (sum(sign_sign)<=sum(sign_sign==F)) {
     Gamma[sign_sign,]=-Gamma[sign_sign,]
     observed_val[subject_ind %in% sl_ind] = - observed_val[subject_ind %in% sl_ind]
     Fitted_values[subject_ind %in% sl_ind] = - Fitted_values[subject_ind %in% sl_ind]
}  else {
sl_ind= l_ind[sign_sign==F]
     Gamma[sign_sign==F,]=-Gamma[sign_sign==F,]
     observed_val[subject_ind %in% sl_ind] = - observed_val[subject_ind %in% sl_ind]
     Fitted_values[subject_ind %in% sl_ind] = - Fitted_values[subject_ind %in% sl_ind]
}

Sign=t(Gamma)
Sign2=cbind(f,Sign)
colnames(Sign2)=c('Freq',as.character(paste('subject',levels(subject))))
Sign2=as.data.frame(Sign2)
save(Sign2,file=paste(path,'Signatures when SeparatePar=',separate,'.RDA'))
rownames(Sign)=as.character(paste('f',f))
colnames(Sign)=as.character(paste('subject',levels(subject)))
write.table(Sign, file=paste(path,"Signatures when SeparatePar=",separate,".txt"), sep="\t", col.names = NA)


colnames(Theta_1)=as.character(paste('subject',levels(subject)))
rownames(Theta_1)=treatment_c
write.table(Theta_1, file=paste(path,"Linear Parameters when SeparatePar=",separate,".txt"), sep="\t", col.names = NA)

f_min=matrix(f_min,ncol=1)
PAR_K=cbind(nlPAR,f_min,1-exp(f_min))
rownames(PAR_K)=rownames(nlPAR)=as.character(paste('subject',levels(subject)))
if (separate==T) {
colnames(nlPAR)= paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)),rep(treatment_c,nPar))
} else {
colnames(nlPAR)=paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)))
}
colnames(PAR_K)=c(c(colnames(nlPAR)),'Func_min','EigValue')
write.table(t(PAR_K), file=paste(path,"Nonlinear Parameters when SeparatePar=",separate,comp_name,".txt"), sep="\t", col.names = NA)

Output_Data=data.frame(subject_ind,Time_ind,treatment_ind, observed_val)

ycoef <- as.data.frame(Gamma)
colnames(ycoef) <- col.name[-c(1:3)]

rownames(Gamma)=levels(subject)
out <- list(call = match.call(),
            ycoef = ycoef, # 'signatures' rownames = subject, colnames = freq
            xcoef = NULL, # yet to implement, should be linear parameters on RHS; rownames - trt levels + cov names + Intercept
            scores = data.frame(subject = as.factor(subject_ind),
                                treatment = as.factor(treatment_ind),# could be NULL?
                                time = Time_ind,
                                xscores = observed_val,
                                yscores = Fitted_values),
            ref = ref, # reference level, not implemented yet (could be NULL?)
            nonlinear.parameters = as.data.frame(nlPAR),
            y = RYr_list,
            x = RFr_list,
            global.roots = PreSmoothingRoots,
            subject.roots = r) # x should be cbind(RFr, G)
class(out) <- "ESLCCA"
out
}

