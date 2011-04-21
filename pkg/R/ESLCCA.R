`ESLCCA` <-
function (DATA,Curve='DoubleExponential',PreSmoothingRoots=0,r=-1,con=0,separate=T,path='c:\\',comp_name,initial.K=NULL) {

# Upload required R libraries
library(MASS)

# Prepare the data
DATA=DATA[!is.na(DATA[,5]),]
col.name=names(DATA)
f=col.name[-c(1:3)]
if (sum(is.na(as.numeric(f)))>0) {f=1:length(f)} else {f=as.numeric(f)}
Y=as.matrix(DATA[,-c(1:3)])
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
Subject=as.factor(as.character(DATA[,names(DATA)=='Subject']))
Treat=as.factor(DATA[,names(DATA)=='Treatment'])
Time= as.numeric(DATA[,names(DATA)=='From'])
nind=nlevels(Subject)
ntreat=nlevels(Treat)
Gamma=matrix(nrow=nind,ncol=length(f))
Fitted_values=Treat_ind=Time_ind= c()
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
        r=min(c(min(table(Subject)),length(f)))
    } else {
        # Default smoothing
        r=c()
        eig=c()
        for (i in 1:nind) {
             Ytemp= scale(Y[Subject==levels(Subject)[i],],scale=F)
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
if (is.null(initial.K)) {
    if (!is.function(Curve)) {
        if (Curve == 'DoubleExponential') {initial.K=rep(c(8.5,9),c(nvar,nvar))}
        if (Curve == 'CriticalExponential') {initial.K=rep(c(8.5),c(nvar))}
    }else{
        stop('No initial values have been specified for the nonlinear parameters.')
    }
}

if (!is.function(Curve)) {
    if (Curve == 'DoubleExponential') {
        Curve <- function(t,Parameters) { K1=Parameters[,1]; K2=Parameters[,2]; exp(-t/exp(K1))-exp(-t/exp(K2)) }
    }
    else {if (Curve == 'CriticalExponential') {
        Curve <- function(t,Parameters) { K1=Parameters[,1]; t*exp(-t/exp(K1)) }
    } }
}

pdf(paste(path,'Double Exponential SeparatePar=',separate, comp_name,'.pdf',sep=""),encoding="ISOLatin2", family="URWHelvetica")
root=1

nlPAR=matrix(nrow=nind,ncol=length(initial.K))

Theta_1=c()
f_min=c()
RYr_list=list()
RFr_list=list()

Subject_ind=c()
corr_mat=matrix(ncol=nind,nrow=nlevels(Treat))
observed_val=c()
for (i in 1:nind) {
     t=c()
     n_vec=c()
     # Do PCA on Y
     IndSubject = Subject==levels(Subject)[i]
     Subject_ind = c(Subject_ind,as.character(Subject)[IndSubject])
     IndV = IndSubject==T & Treat=='Control'
     nv = sum(IndV); tv = Time[IndV]
     Yr = Y[IndV,]
     Treat_c=as.character(levels(Treat))
     Treat_c=Treat_c[Treat_c!='Control']
     char.Treat=as.character(Treat)
     Treat_ind= c(Treat_ind,char.Treat[IndV])
     F0=rep(rep(0,length(Treat_c)),nv)
     for (j in 1:length(Treat_c))  {
     IndL = IndSubject==T & Treat==Treat_c[j]
     nl = sum(IndL); tl = Time[IndL]
     Yr = rbind(Yr,Y[IndL,])
     n_vec=c(n_vec,nl)
     t=c(t,tl)
     Treat_ind= c(Treat_ind,char.Treat[IndL])
     new_F=rep(0,length(Treat_c))
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
     obj.f=function(Kvector) {
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

l_ind=levels(Subject)
mydata.pca <- princomp(t(Gamma))
loading <- mydata.pca$loadings[,1]
sign_sign= loading>0
sl_ind=l_ind[sign_sign]
if (sum(sign_sign)<=sum(sign_sign==F)) {
     Gamma[sign_sign,]=-Gamma[sign_sign,]
     observed_val[Subject_ind %in% sl_ind] = - observed_val[Subject_ind %in% sl_ind]
     Fitted_values[Subject_ind %in% sl_ind] = - Fitted_values[Subject_ind %in% sl_ind]
}  else {
sl_ind= l_ind[sign_sign==F]
     Gamma[sign_sign==F,]=-Gamma[sign_sign==F,]
     observed_val[Subject_ind %in% sl_ind] = - observed_val[Subject_ind %in% sl_ind]
     Fitted_values[Subject_ind %in% sl_ind] = - Fitted_values[Subject_ind %in% sl_ind]
}

Sign=t(Gamma)
Sign2=cbind(f,Sign)
colnames(Sign2)=c('Freq',as.character(paste('Subject',levels(Subject))))
Sign2=as.data.frame(Sign2)
save(Sign2,file=paste(path,'Signatures when SeparatePar=',separate,'.RDA'))
rownames(Sign)=as.character(paste('f',f))
colnames(Sign)=as.character(paste('Subject',levels(Subject)))
write.table(Sign, file=paste(path,"Signatures when SeparatePar=",separate,".txt"), sep="\t", col.names = NA)


colnames(Theta_1)=as.character(paste('Subject',levels(Subject)))
rownames(Theta_1)=Treat_c
write.table(Theta_1, file=paste(path,"Linear Parameters when SeparatePar=",separate,".txt"), sep="\t", col.names = NA)

f_min=matrix(f_min,ncol=1)
PAR_K=cbind(nlPAR,f_min,1-exp(f_min))
rownames(PAR_K)=rownames(nlPAR)=as.character(paste('Subject',levels(Subject)))
if (separate==T) {
colnames(nlPAR)= paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)),rep(Treat_c,nPar))
} else {
colnames(nlPAR)=paste(rep(paste('Parameter',1:nPar),rep(nvar,nPar)))
}
colnames(PAR_K)=c(c(colnames(nlPAR)),'Func_min','EigValue')
write.table(t(PAR_K), file=paste(path,"Nonlinear Parameters when SeparatePar=",separate,comp_name,".txt"), sep="\t", col.names = NA)

Output_Data=data.frame(Subject_ind,Time_ind,Treat_ind, observed_val)

dev.off()
rownames(Gamma)=levels(Subject)
out <- list(frequencies=f,signatures=Gamma, fitted.values = observed_val,
            y = Fitted_values, # observed_val & Fitted_values currently defined opposite way from what you'd expect!
            subject = as.factor(Subject_ind), time = Time_ind, treatment = as.factor(Treat_ind),
            NonLinearPars=nlPAR,RYr=RYr_list,RFr=RFr_list)
class(out) <- "ESLCCA"
out
}

