##
##fn_tsbugs1.R: uv models
##fn_tsbugs2.R: general functions
##fn_tsbugs4.R: garch models
##fn_tsbugs4.R: mv models
##

##
##GARCH
##
#y=svpdx$pdx; ar.order=1; sim=TRUE; k=10; beg=2; mean.centre=TRUE; beg=ar.order+1; ar.prior="dnorm(0,1)"; mean.prior="dnorm(0,1)";
#arch.order=1; garch.order=1; vol.beg=beg+max(arch.order,garch.order); 
#vol.mean.prior="dgamma(0.000001,0.000001)"; arch.prior="dgamma(0.000001,0.000001)"; garch.prior="dgamma(0.000001,0.000001)"; space=FALSE
garch.bugs<-function(y, ar.order=0, k=NULL, sim=FALSE, 
                     mean.centre=FALSE, beg=ar.order+1,
                     mean.prior=ar.prior, ar.prior="dnorm(0,1)",
                     arch.order=1, garch.order=1, vol.beg=max(arch.order,garch.order)+1,
                     h0.prior="dgamma(0.000001,0.000001)", 
                     vol.mean.prior="dgamma(0.000001,0.000001)", 
                     arch.prior="dgamma(0.000001,0.000001)", 
                     garch.prior="dgamma(0.000001,0.000001)", 
                     space=FALSE){
  y<-c(y)
  n<-length(y)
  if(!is.null(k)){
    y<-c(y,rep(NA,k))
  }
  k<-length(y)-max(which(!is.na(y)))
  if(beg<ar.order)
    stop("The value of beg must be at least 1 greater than the number of lags")
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",beg,":",n+k,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2[t])",
         "\tisigma2[t] <- 1/h[t]",
         "}")
  bug<-c(bug, lik)
  #ymean
  ymean<-c("#mean",
           paste0("for(t in ",beg,":",n+k,"){"),
           y.mean<-c("\ty.mean[t] <- 0",
                     "}")
  )
  if(ar.order==0 & mean.centre==TRUE)  ymean[3]<-"\ty.mean[t] <- phi0"
  if(ar.order!=0 & mean.centre==FALSE)  ymean[3]<-paste0("\ty.mean[t] <- ",paste0("phi",1:ar.order,"*y[t-",1:ar.order,"]",collapse=" + "))
  if(ar.order!=0 & mean.centre==TRUE)  ymean[3]<-paste0("\ty.mean[t] <- phi0 + ",paste0("phi",1:ar.order,"*(y[t-",1:ar.order,"]-phi0)",collapse=" + "))
  bug<-c(bug, ymean)
  
  #hmean
  hmean<-c("#volatility",
           paste0("for(t in ",beg,":",vol.beg-1,"){"),
           paste0("\th[t] ~ ",h0.prior),
           "}",
           paste0("for(t in ",vol.beg,":",n+k,"){"),
           paste0("\th[t] <- psi0 + ",paste0("omega",1:arch.order,"*pow(y[t-",1:arch.order,"],2)",collapse=" + "),
                             " + ",paste0("psi",1:garch.order,"*(h[t-",1:garch.order,"])",collapse=" + ")),
           "}",
           "")
  bug<-c(bug, hmean)
  
  #priors
  ar.priors<-paste0("phi",0:ar.order," ~ ",ar.prior)
  ar.priors<-ar.priors[-1]
  if(mean.centre==TRUE)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
  psi0.priors<-paste0("psi0 ~ ",vol.mean.prior)
  arch.priors<-paste0("omega",1:arch.order," ~ ",arch.prior)
  garch.priors<-paste0("psi",1:garch.order," ~ ",garch.prior)
  vol.priors<-c(psi0.priors,arch.priors, garch.priors)
  bug<-c(bug,"#priors",ar.priors,vol.priors,"")
  
  #forecast
  forc<-NULL
  if(k!=0){
    forc<-c("#forecast",
            paste("for(t in ",n+1,":",n+k,"){",sep=""),
            "\ty.new[t] <- y[t]",
            "}",
            "")
    bug<-c(bug,forc)
  }
  
  #simulation
  if(sim==TRUE){
    ysim<-c("#simulation",
            paste("for(t in ",beg,":",n,"){",sep=""),
            "\ty.mean.c[t] <- cut(y.mean[t])",
            "\tisigma2.c[t] <- cut(isigma2[t])",
            "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c[t])",
            "}",
            "")
    bug<-c(bug,ysim)
  }
  bug<-c(bug,"}","")
  
  if(space==FALSE){
    bug<-bug[-nchar(bug)!=0]
    if(length(grep("#mean", bug))>0)
      bug<-bug[-grep("#mean", bug)]
    if(length(grep("#volatility", bug))>0)
      bug<-bug[-grep("#volatility", bug)]
  }
  
  p1<-grep("#likelihood",bug)
  p2<-grep("#prior",bug)
  if(k!=0 & sim==TRUE){
    p3<-grep("#forecast",bug); p4<-grep("#simulation",bug)
  }
  if(k!=0 & sim==FALSE){
    p3<-grep("#forecast",bug); p4<-length(bug)
  }
  if(k==0 & sim==TRUE){
    p3<-grep("#simulation",bug); p4<-p3
  } 
  if(k==0 & sim==FALSE){
    p3<-length(bug); p4<-p3
  } 
  p5<-length(bug)
  
  bug<-list(bug=bug,
            data=list(y=y),
            info=list(n=n,k=k,nh=n+k,
                      args=mget(names(formals()),sys.frame(sys.nframe()))[-1],
                      variance="GARCH",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecast=NULL,
                      simulation=NULL))
  if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
}

#garch.bugs(y,k=5, ar.order=4,sim=TRUE, arch.order=3, mean.centre=T, beg=10)


