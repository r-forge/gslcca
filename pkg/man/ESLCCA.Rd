\name{ESLCCA}
\alias{ESLCCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ESLCCA }
\description{
  Detecting drug entry into the brain is important clinically.  This functions detects and apply Semi Linear Canonical Correlation (SLCCA) to provide the PK/PD-like signals in EEG data.  The algorithm can be applied in a regression-like framework where a range of different models can be readily fitted to obtain drug “signatures” and smooth time courses.  This function uses a data smoothing approach based on Principal Components Analysis (PCA) to obtain smoothed drug “signatures” which can be readily compared between different subjects.
}
\usage{
ESLCCA(DATA, Curve = "DoubleExponential", PreSmoothing = 0, r = -1, con = 0, separate = T, path = "c:\\", comp_name, initial.K = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{ data frame with the first column called 'Subject' (indicates the subject's ID), the second called 'Treatment' (indicates the type of treatment for each subject), the third called 'From' (contains the lower limit of the time interval for the FFT) and the rest columns contain the power spectrum at each frequency. So the number of columns is 3 plus the number of frequencies. }
  \item{Curve}{ either an R function with the curve to be fitted or a character string. If a character string then it can be one of the two, 'DoubleExponential' (the default) or 'CriticalExponential'. }
  \item{PreSmoothingRoots}{ a scalar. If it is greater than 0 then smoothing is performed using PCA and the number of roots is equal to PreSmoothingRoots. PreSmoothingRoots can also take the value -1 if we want the algorithm itself to determine the number of components by finding such k that l_k/l_1 >= 10.E-3 and S_{j=1}^k l_j/S_{j=1}^r l_j >= 0.98, where k=1,...,r, r=rank(DATA[,-c(1:3)]) and l_1, l_2,... are the eigenvalues of the variance matrix. Finally if 'roots' takes the number 0 (the default) then no pre-smoothing is performed. Pre-smoothing the data matrix (for all subjects silmutaneously) is a way to remove artifacts.  }
  \item{r}{ ~~Describe \code{r} here~~ }
  \item{con}{ takes the values 0 or 1. If 1 then no ajustment for the overall mean is performed, i.e R=diag(I). If 0 then we ajust for the mean which is the default, i.e . R = diag(I)-I%*%solve(t(I)%*%I)%*%t(I) }
  \item{separate}{ if TRUE then the parameters (linear and non-linear) are allowed to be different for different groups of observations. For example if we have EEG data for a subject recorded at different periods then we can assume different characteristics for the different periods. }
  \item{path}{ the directory where to save the results (all generated plots in pdf and all the results in text files and rda. ) }
  \item{comp_name}{ a character string with the name of the compound (it is optional, it appears at the name of the saved files and its purpose is to help discriminate among saved files) }
  \item{initial.K}{ a vector with the initial values for each of the nonlinear parameters. It is optional only when Curve='DoubleExponential' or Curve='CriticalExponential'. However it is preferred for the user to specify the initial values since good initial values are essential and are always case dependant.  }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (DATA,Curve='DoubleExponential',PreSmoothing=list(Smoothing=F,roots=0),r=-1,con=0,separate=T,path='c:\\',comp_name,initial.K=NULL) {

# Upload required R libraries
library(MASS)

# Prepare the data
DATA=DATA[!is.na(DATA[,5]),]
col.name=names(DATA)
f=col.name[-c(1:3)]
if (sum(is.na(as.numeric(f)))>0) {f=1:length(f)} else {f=as.numeric(f)}
Y=as.matrix(DATA[,-c(1:3)])
# Pre-smoothing of the complete data set for artifacts removal
if (PreSmoothing$Smoothing==T) {
    if  (PreSmoothing$roots==-1) {
         Ytemp= scale(Y,scale=F)
         eig=cbind(eigen(t(Ytemp)\%*\%Ytemp)$values)
         nroots=1:length(eig)
         PreSmoothing$roots=max(nroots[(eig/eig[1])>0.001 & cumsum(eig)/sum(eig)<0.98]) + 1
    }
    SVD = svd( Y, nu=PreSmoothing$roots, nv=PreSmoothing$roots )
    Y = SVD$u\%*\%diag(SVD$d[1:PreSmoothing$roots])\%*\%t(SVD$v)
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
             eig=cbind(eigen(t(Ytemp)\%*\%Ytemp)$values)
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
     Yr = SVD$u\%*\%diag(SVD$d[1:r])\%*\%t(SVD$v)
     # Calculate the idempotent matrix R
     I = rep(1,nr)
     if (con==0) {
         R = diag(I)-I\%*\%solve(t(I)\%*\%I)\%*\%t(I)
     } else {
         R = diag(I)
     }
     RYr=R\%*\%Yr
     # Reduce the dimensionality using svd
     rank_RFr=qr(RYr)$rank
     SVD = svd( RYr, nu=rank_RFr, nv=rank_RFr )
     D<-diag(SVD$d[1:rank_RFr])
     U1<-SVD$u
     V<-SVD$v
     RYr=U1\%*\%D
     # Define F
     F0=matrix(F0,nrow=nr,byrow=T)
     Fnl=function(Parameters) {
         Fb=F0
         Fb[F0==1]=Curve(t,Parameters)
         Fb
         }
     S11=t(RYr)\%*\%RYr
     S11.eig <- eigen(S11)
     S11.val=suppressWarnings(sqrt(S11.eig$values))
     S11.val[is.na(S11.val)]=0
     S11.sqrt <- ginv(S11.eig$vectors \%*\% diag(S11.val) \%*\% t(S11.eig$vectors))
     nPar=length(initial.K)/length(size())  # computes the number of parameters in 'Curve'
     temp1=rep(size(),nPar)
     obj.f=function(Kvector) {
           RFr = R\%*\%Fnl(matrix(rep(Kvector,temp1),ncol=nPar))
           S12=t(RYr)\%*\%RFr
           log(1-eigen(S11.sqrt\%*\%S12\%*\%ginv(t(RFr)\%*\%RFr)\%*\%t(S12)\%*\%S11.sqrt)$values[1])
     }
     # Minimize ln(1-largesteigenvalue(Corr))
     nlPAR[i,]=optim(initial.K,obj.f,control=list(warnOnly = T))$par
     RFr = R\%*\%Fnl(matrix(rep(nlPAR[i,],temp1),ncol=nPar))
     S12=t(RYr)\%*\%RFr
     S22.inv.S21=ginv(t(RFr)\%*\%RFr)\%*\%t(S12)
     Eigvectors=eigen(S11.sqrt\%*\%S12\%*\%S22.inv.S21\%*\%S11.sqrt)$vectors
     optEigvalue=Eigvectors[,root]
     Gamma[i,]= V\%*\%S11.sqrt\%*\%optEigvalue
     Gamma= Re(Gamma)
     #Theta= S22.inv.S21\%*\%Gamma[i,]
     Theta= Re(S22.inv.S21\%*\%Eigvectors)
     RFrTheta=RFr\%*\%Theta
     Norm=apply(RFrTheta^2,2,sum)
     Theta=Theta/sqrt(matrix(Norm,nrow(Theta),ncol=ncol(Theta),byrow=T))
     RFrTheta=RFr\%*\%Theta[,root]
     hRFrTheta=RFrTheta[-c(1:(length(t)+10))]
     cond=hRFrTheta[abs(hRFrTheta)==max(abs(hRFrTheta))]
     if   (cond[1]<0) {Gamma[i,]=-Gamma[i,]; RFrTheta=-RFrTheta}
     #if   (min(abs(RFrTheta))>max(abs(RFrTheta))) {Gamma[i,]=-Gamma[i,]; RFrTheta=-RFrTheta}
     par(mfrow=c(2,1),cex=0.8)
     plot(f,Gamma[i,],xlab='Frequency',ylab='Signature',type='l')
     title(paste('Signature for subject',levels(Subject)[i]))
     Fit_values=RYr\%*\%t(V)\%*\%Gamma[i,]
     Fitted_values=c(Fitted_values,Fit_values)
     observed_val=c(observed_val,RFrTheta)
     plot(c(tv,t),RFrTheta,ylim=range(c(RFrTheta,Fit_values)),pch='*',xlab='Time',ylab='Fitted values',col=rep(1:(length(Treat_c)+1),c(nv,n_vec)))
     points(c(tv,t),Fit_values,pch='.',col=rep(1:(length(Treat_c)+1),c(nv,n_vec)))
     title(paste('Fitted values for subject',levels(Subject)[i]))
     legend(x="topright",y=NULL, c('Control',Treat_c), col =1:(length(Treat_c)+1), text.col = "green4", lty = 1, pch = NULL, merge = TRUE)

     Theta_1=cbind(Theta_1,Theta[,1])
     f_min=c(f_min,Re(obj.f(nlPAR[i,])))
     RYr_list[[i]]=RYr\%*\%t(V)
     RFr_list[[i]]=RFr
  }

l_ind=levels(Subject)
mydata.pca <- princomp(t(Gamma))
loading <- mydata.pca$loadings[,1]
sign_sign= loading>0
sl_ind=l_ind[sign_sign]
if (sum(sign_sign)<=sum(sign_sign==F)) {
     Gamma[sign_sign,]=-Gamma[sign_sign,]
     observed_val[Subject_ind \%in\% sl_ind] = - observed_val[Subject_ind \%in\% sl_ind]
     Fitted_values[Subject_ind \%in\% sl_ind] = - Fitted_values[Subject_ind \%in\% sl_ind]
  }  else {
sl_ind= l_ind[sign_sign==F]
     Gamma[sign_sign==F,]=-Gamma[sign_sign==F,]
     observed_val[Subject_ind \%in\% sl_ind] = - observed_val[Subject_ind \%in\% sl_ind]
     Fitted_values[Subject_ind \%in\% sl_ind] = - Fitted_values[Subject_ind \%in\% sl_ind]
  }

par(mfrow=c(1,1),cex=0.8)
matplot(f,t(Gamma),col=1:nind,type='l',xlab='Frequency',ylab='Signature')
title('Signatures corresponding to different subjects')
legend(x="topright",y=NULL, paste('Subject',levels(Subject)), col =1:nind, text.col = "green4", lty = 1:nind, pch = NULL, merge = TRUE)

plot(apply(Gamma,2,mean),xlab='Frequency',type='l')
title('Mean Signature' )


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

Output_Data=data.frame(Subject_ind,Time_ind,Treat_ind,Fitted_values,observed_val)

# Calculate the mean fitted curves
data=Output_Data[c(1:3,5)]
time_treat=paste(data[,2],data[,3])
Subject_ind=as.numeric(data$Subject_ind)
observed_val=as.numeric(data$observed_val)
data=data.frame(time_treat,Subject_ind,observed_val)
wide <- reshape(data, v.names="observed_val", idvar="time_treat",timevar="Subject_ind", direction="wide")
D=matrix(as.numeric(as.matrix(wide)),nrow=nrow(wide),ncol=ncol(wide))
Time=as.numeric(levels(as.factor(Output_Data$Time_ind)))
# Plot the mean fitted curves
plot(rep(Time,length(Treat_c)+1),apply(D,1,mean,na.rm=TRUE),col=rep(1:(length(Treat_c)+1),rep(length(Time),length(Treat_c)+1)),xlab='Time',ylab='Fitted values')
title('Mean fitted curves')
legend(x="topright",y=NULL, c('Control',Treat_c), col =1:(length(Treat_c)+1), text.col = "green4", lty = 1, pch = NULL, merge = TRUE)


dev.off()
rownames(Gamma)=levels(Subject)
list(frequencies=col.name[-c(1:3)],Signature=Gamma,fit_obs=Output_Data ,NonLinearPars=nlPAR,RYr=RYr_list,RFr=RFr_list)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
