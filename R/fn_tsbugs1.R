##
##fn_tsbugs1.R: uv models
##fn_tsbugs2.R: general functions
##fn_tsbugs4.R: garch models
##fn_tsbugs4.R: mv models
##

##
##AR Model
##
#y=diff(lh); ar.order=1; sim=TRUE; k=NULL; beg=ar.order+1; mean.centre=FALSE; beg=ar.order+1; ar.prior="dnorm(0,1)"; tol.prior="dgamma(0.000001,0.000001)";
#ar.prior="dnorm(0,1)"; tol.prior="dgamma(0.000001,0.000001)"; var.prior=NULL; sd.prior=NULL; mean.prior=ar.prior
ar.bugs<-function(y, ar.order=1, k=NULL, sim=FALSE, 
                  mean.centre=FALSE, beg=ar.order+1,
                  mean.prior=ar.prior, ar.prior="dnorm(0,1)", tol.prior="dgamma(0.000001,0.000001)", var.prior=NULL, sd.prior=NULL,
                  space=FALSE){
  y<-c(y)
  n<-length(y)
  if(!is.null(k)){
    y<-c(y,rep(NA,k))
  }
  k<-length(y)-max(which(!is.na(y)))
  if(beg<ar.order)
    stop("The value of beg must be at least 1 greater than the number of lags")
  if(!is.null(var.prior) | !is.null(sd.prior)){
    tol.prior<-NULL
  }
  if(length(c(tol.prior,var.prior,sd.prior))>1)
    stop("Only one of tol.prior, var.prior or sd.prior should be given. Set others to null")
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",beg,":",n+k,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2)",
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
  bug<-c(bug, ymean, "")
  
  #prior
  if(!is.null(tol.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("isigma2 ~ ",tol.prior),
                 "sigma <- pow(isigma2,-0.5)")
  }
  if(!is.null(var.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma2 ~ ",var.prior),
                 "isigma2 <- pow(sigma2,-1)")
  }
  if(!is.null(sd.prior)){
    ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
                 paste0("sigma ~ ",sd.prior),
                 "isigma2 <- pow(sigma,-2)")
  }
  ar.priors<-ar.priors[-1]
  if(mean.centre==TRUE)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
  bug<-c(bug,"#priors",ar.priors,"")
  
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
  ysim<-NULL
  if(sim==TRUE){
    ysim<-c("#simulation",
            "isigma2.c <- cut(isigma2)",
            paste("for(t in ",beg,":",n,"){",sep=""),
            "\ty.mean.c[t] <- cut(y.mean[t])",
            "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c)",
            "}",
            "")
    bug<-c(bug,ysim)
  }
  bug<-c(bug,"}","")
  #print.tsbugs(list(bug=bug))
  
  if(space==FALSE){
    bug<-bug[-nchar(bug)!=0]
    if(length(grep("#mean", bug))>0)
      bug<-bug[-grep("#mean", bug)]
  }
    
  p1<-grep("#likelihood",bug)
  p2<-grep("#priors",bug)
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
                      variance="CV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecast=NULL,
                      simulation=NULL))
  if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
}



##
##SV
##
#y=diff(lh); ar.order=1; sim=TRUE; k=10; beg=2; mean.centre=TRUE; beg=ar.order+1; ar.prior="dnorm(0,1)"; mean.prior="dnorm(0,1)";
#sv.order=1; sv.mean.prior1="dgamma(0.000001,0.000001)"; sv.ar.prior1="dbeta(1,1)"; sv.tol.prior="dgamma(0.01,0.01)";  sv.mean.prior2=NULL; sv.ar.prior2=NULL
sv.bugs<-function(y, ar.order=0, k=NULL, sim=FALSE, 
                  mean.centre=FALSE, beg=ar.order+1,
                  mean.prior=ar.prior, ar.prior="dnorm(0,1)",
                  sv.order=1, sv.beg=beg+sv.order,
                  sv.mean.prior1="dnorm(0,0.001)", sv.mean.prior2=NULL,
                  sv.ar.prior1="dunif(0,1)", sv.ar.prior2=NULL,
                  sv.tol.prior="dgamma(0.01,0.01)",
                  space=FALSE){
  y<-c(y)
  n<-length(y)
  if(!is.null(k)){
    y<-c(y,rep(NA,k))
  }
  k<-length(y)-max(which(!is.na(y)))
  if(beg<ar.order)
    stop("The value of beg must be at least 1 greater than the number of lags")
  if(!is.null(sv.ar.prior2)){
    sv.ar.prior1<-NULL
  }
  if(!is.null(sv.mean.prior2)){
    sv.mean.prior1<-NULL
  }
  if(length(c(sv.mean.prior1,sv.mean.prior2))>1)
    stop("Only one of sv.mean.prior1 or sv.mean.prior2 should be given. Set others to null")
  if(length(c(sv.ar.prior1,sv.ar.prior2))>1)
    stop("Only one of sv.ar.prior1 or sv.ar.prior2 should be given. Set others to null")  
  
  bug<-c("model{","")
  #likelihood
  lik<-c("#likelihood",
         paste0("for(t in ",beg,":",n+k,"){"),
         "\ty[t] ~ dnorm(y.mean[t], isigma2[t])",
         "\tisigma2[t] <- exp(-h[t])",
         "\th[t] ~ dnorm(h.mean[t], itau2)",
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
           paste0("for(t in ",beg,":",sv.beg-1,"){"),
           "\th.mean[t] <- psi0",
           "}",
           paste0("for(t in ",sv.beg,":",n+k,"){"),
           paste0("\th.mean[t] <- psi0 + ",paste0("psi",1:sv.order,"*(h[t-",1:sv.order,"]-psi0)",collapse=" + ")),
           "}",
           "")
  bug<-c(bug, hmean)
  
  #priors
  ar.priors<-paste0("phi",0:ar.order," ~ ",ar.prior)
  ar.priors<-ar.priors[-1]
  if(mean.centre==TRUE)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
  if(!is.null(sv.mean.prior2)){
    sv.priors<-c(paste0("psi0.star ~ ",sv.mean.prior2),
                 "psi0 <- -log(psi0.star)")
  }
  if(!is.null(sv.mean.prior1)){
    sv.priors<-paste0("psi0 ~ ",sv.mean.prior1)
  }
  if(!is.null(sv.ar.prior2)){
    sv.priors<-c(sv.priors,
                 paste0("psi",1:sv.order," ~ ",sv.ar.prior2))
  }
  if(!is.null(sv.ar.prior1)){
    sv.priors<-c(sv.priors,
                 paste0("psi",1:sv.order,".star ~ ",sv.ar.prior1),
                 paste0("psi",1:sv.order," <- 2*psi",1:sv.order,".star-1"))
  }
  sv.priors<-c(sv.priors, 
               paste0("itau2 ~ ",sv.tol.prior),
               "tau <- pow(itau2,-0.5)")
  bug<-c(bug,"#priors",ar.priors,sv.priors,"")
  
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
                      variance="SV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecast=NULL,
                      simulation=NULL))
  if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
}
#sv.bugs(y,k=5, ar.order=4,sim=T, sv.order=10, mean.centre=T, beg=10)


##
##RV
##
#y=diff(lh); ar.order=1; sim=TRUE; k=10; beg=2; mean.centre=TRUE; beg=ar.order+1; ar.prior="dnorm(0,1)"; 
# rv.tol0.prior="dgamma(0.000001,0.000001)"; rv.eps.prior="dbeta(1, 100)"; rv.var.prior="dgamma(0.01,0.01)"
rv.bugs<-function(y, ar.order=0, k=NULL, sim=FALSE, 
                  mean.centre=FALSE, beg=ar.order+1,
                  mean.prior=ar.prior, 
                  ar.prior="dnorm(0,1)",
                  rv.tol0.prior="dgamma(0.000001,0.000001)", 
                  rv.eps.prior="dbeta(1, 100)", 
                  rv.ilambda2.prior="dgamma(0.01,0.01)",
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
         "\tisigma2[t] <- exp(-h[t])",
         "\th[t] <- 2*lsig[t]",
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
  
  #rv
  rv<-c("#variance",
        paste0("lsig[",beg,"] <- -0.5*log(isig02)"),
        paste0("for(t in ",beg+1,":",n+k,"){"),
        "\tlsig[t] <- lsig[t-1]+(delta[t]*beta[t])",
        "\tdelta[t] ~ dbern(epsilon)",
        "\tbeta[t] ~ dnorm(0,ilambda2)",
        "}",
        "")
  bug<-c(bug, rv)
  
  #priors
  ar.priors<-paste0("phi",0:ar.order," ~ ",ar.prior)
  ar.priors<-ar.priors[-1]
  if(mean.centre==TRUE)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)  
  rv.priors<-c(paste0("isig02 ~ ",rv.tol0.prior),
               "sig0 <- pow(isig02,-0.5)",
               paste0("epsilon ~ ",rv.eps.prior),
               paste0("ilambda2 ~ ",rv.ilambda2.prior),
               "lambda <- pow(ilambda2,-0.5)")
  bug<-c(bug,"#priors",ar.priors,rv.priors,"")
  
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
    if(length(grep("#variance", bug))>0)
      bug<-bug[-grep("#variance", bug)]
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
                      variance="RV",
                      likelihood=p1:(p2-1),
                      priors=p2:(p3-1),
                      forecast=NULL,
                      simulation=NULL))
  if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
  if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
  class(bug)<-"tsbugs"
  return(bug)
  class(bug)<-"tsbugs"
  return(bug)
}



##
##random walk (with drift)
##
# rw.bugs<-function(n, k=0, simulations=FALSE, drift=FALSE, beg=1,
#                   drift.prior="dnorm(0,1)", tol.prior="dgamma(0.000001,0.000001)"){
#   bug<-c("model{","")
#   #likelihood
#   lik<-c("#likelihood",
#          paste("for(t in ",beg,":",n+k,"){",sep=""),
#          "\ty[t] ~ dnorm(y.mean[t], isigma2)",
#          "}",
#          paste("for(t in ",beg+1,":",n+k,"){",sep=""))
#   bug<-c(bug, lik)
#   if(drift==FALSE)  y.mean<-c("\ty.mean[t] <- y[t-1]","}")
#   if(drift==TRUE)  y.mean<-c("\ty.mean[t] <- mu + y[t-1]","}")
#   bug<-c(bug, y.mean, "")
#   #priors
#   rw.priors<-c(paste("isigma2 ~ ",tol.prior,sep=""),
#                "sigma <- pow(isigma2,-0.5)")
#   if(drift==FALSE) bug<-c(bug,"#priors",rw.priors,"")
#   if(drift==TRUE) bug<-c(bug,"#priors",paste("mu ~ ",drift.prior,sep=""),rw.priors,"")
#   #forecast
#   if(k!=0){
#     forc<-c("#forecast",
#             paste("for(t in ",n+1,":",n+k,"){",sep=""),
#             "\ty.new[t] <- y[t]",
#             "}",
#             "")
#     bug<-c(bug,forc)
#   }  
#   #simulation
#   if(simulations==TRUE){
#     ysim<-c("#simulation",
#             "isigma2.c <- cut(isigma2)",
#             paste("for(t in ",beg,":",n,"){",sep=""),
#             "\ty.mean.c[t] <- cut(y.mean[t])",
#             "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c)",
#             "}",
#             "")
#     bug<-c(bug,ysim)
#   }
#   bug<-c(bug,"}","")
#   class(bug)<-"mbugs"
#   bug
# }
# rw.bugs(n=20,k=5, sim=T, drift=T)



##
##MA Model
##
# ma.bugs<-function(y, ma.order=1, n=length(y)-k, k=length(y)-max(which(!is.na(y))), b=min(which(!is.na(y)))-1, simulations=FALSE, 
#                   mean.centre=FALSE, beg=ma.order+1,
#                   mean.prior=ma.prior, ma.prior="dnorm(0,1)", tol.prior="dgamma(0.000001,0.000001)", var.prior=NULL, sd.prior=NULL){
#   bug<-c("model{","")
#   #likelihood
#   lik<-c("#likelihood",
#          paste("for(t in ",max(1,b-beg),":",b+n+k,"){",sep=""),
#          "\ty[t] ~ dnorm(y.mean[t], isigma2)",
#          "\te[t] ~ y[t] - y.mean[t]",
#          "}")
#   if(b!=0 & b>beg)  lik[2]<-paste("for(t in 1:",b+n+k,"){",sep="")
#   bug<-c(bug, lik)
#   #latent data mean
#   int<-c("#priors for latent data",
#          paste("for(t in 1:",max(beg,b),"){",sep=""),
#          y.mean<-c("\ty.mean[t] <- 0",
#                    "}")
#   )
#   if(mean.centre==TRUE)  int[3] <- c("\ty.mean[t] <- phi0")
#   if(b!=0 & b>beg)  bug<-c(bug, int)
#   #ymean
#   ymean<-c("#data mean",
#            paste("for(t in ",max(beg,b+1),":",b+n+k,"){",sep=""),
#            y.mean<-c("\ty.mean[t] <- 0",
#                      "}")
#   )
#   if(ar.order==0 & mean.centre==TRUE)  ymean[3]<-"\ty.mean[t] <- phi0"
#   if(ar.order!=0 & mean.centre==FALSE)  ymean[3]<-paste("\ty.mean[t] <- ",paste("phi",1:ar.order,"*y[t-",1:ar.order,"]",sep="",collapse=" + "),sep="")
#   if(ar.order!=0 & mean.centre==TRUE)  ymean[3]<-paste("\ty.mean[t] <- phi0 + ",paste("phi",1:ar.order,"*(y[t-",1:ar.order,"]-phi0)",sep="",collapse=" + "),sep="")
#   bug<-c(bug, ymean, "")
#   #priors
#   if(!is.null(tol.prior)){
#     ar.priors<-c(paste("phi",1:ar.order," ~ ",ar.prior,sep=""), 
#                  paste("isigma2 ~ ",tol.prior,sep=""),
#                  "sigma <- pow(isigma2,-0.5)")
#   }
#   if(!is.null(var.prior)){
#     ar.priors<-c(paste("phi",1:ar.order," ~ ",ar.prior,sep=""), 
#                  paste("sigma2 ~ ",var.prior,sep=""),
#                  "isigma2 <- pow(sigma2,-1)")
#   }
#   if(!is.null(sd.prior)){
#     ar.priors<-c(paste("phi",1:ar.order," ~ ",ar.prior,sep=""), 
#                  paste("sigma ~ ",sd.prior,sep=""),
#                  "isigma2 <- pow(sigma,-2)")
#   }
#   if(mean.centre==FALSE) bug<-c(bug,"#priors",ar.priors,"")
#   if(mean.centre==TRUE) bug<-c(bug,"#priors",paste("phi0 ~ ",mean.prior,sep=""),ar.priors,"")
#   #backcast
#   if(b!=0 & b>beg){
#     back<-c("#backcasts",
#             paste("for(t in 1:",b,"){",sep=""),
#             "\ty.old[t] <- y[t]",
#             "}",
#             "")
#     bug<-c(bug,back)
#   }
#   #forecast
#   if(k!=0){
#     forc<-c("#forecast",
#             paste("for(t in ",n+1,":",n+k,"){",sep=""),
#             "\ty.new[t] <- y[t]",
#             "}",
#             "")
#     bug<-c(bug,forc)
#   }  
#   #simulation
#   if(simulations==TRUE){
#     ysim<-c("#simulation",
#             "isigma2.c <- cut(isigma2)",
#             paste("for(t in ",max(beg,b+1),":",n,"){",sep=""),
#             "\ty.mean.c[t] <- cut(y.mean[t])",
#             "\ty.sim[t] ~ dnorm(y.mean.c[t],isigma2.c)",
#             "}",
#             "")
#     bug<-c(bug,ysim)
#   }
#   bug<-c(bug,"}","")
#   class(bug)<-"mbugs"
#   bug
# }
# 

# m1<-ar.bugs(y, ar.order=4,b=10)
# print(m1)
# bug<-m1
# nodes<-function(bug=NULL){
#   if(class(bug)!="tsbugs")
#     stop("bug must be of class bug (i.e. created using function in ts4BUGS")
#   forc<-NULL;back<-NULL
#   lik<-bug[(grep("#likelihood",bug)+2):(grep("^$",bug)[grep("^$",bug)>grep("#likelihood",bug)][1]-2)]
#   prior<-bug[(grep("#priors",bug)+1):(grep("^$",bug)[grep("^$",bug)>grep("#priors",bug)][1]-1)]
#   if(length(grep("#backcasts",bug))>0)
#     back<-bug[(grep("#backcasts",bug)+2):(grep("^$",bug)[grep("^$",bug)>grep("#backcasts",bug)][1]-2)]
#   if(length(grep("#forcasts",bug,value=0))>0)
#     forc<-bug[(grep("#forecast",bug)+2):(grep("^$",bug)[grep("^$",bug)>grep("#forecast",bug)][1]-2)]
#   
#   nds<-ar.priors
#   
#   if(!is.null(forc))
#     forc<-clean(forc)
#   if(!is.null(back))
#     back<-clean(back)
#   
#   list(prior=clean(prior),forc=forc,back=back,lik=clean(lik))
# }
# 
