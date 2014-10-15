# ##
# ##fn_tsbugs1.R: uv models
# ##fn_tsbugs2.R: general functions
# ##fn_tsbugs4.R: garch models
# ##fn_tsbugs4.R: mv models
# ##
# 
# anzdx<-read.csv("AUDNZD.csv", header=FALSE)
# names(anzdx)<-c("adx","nzdx")
# anzdx<-ts(anzdx,start=c(1994,1), frequency=52)
# plot(anzdx, mar.multi = c(0,0,0,0), oma.multi = c(2,2,0.1,0.1))
# 
# 
# solve(cov(anzdx))
# 
# y<-anzdx
# ##
# ##AR Model
# ##
# #y=diff(lh); ar.order=1; sim=TRUE; k=NULL; beg=ar.order+1; mean.centre=FALSE; beg=ar.order+1; ar.prior="dnorm(0,1)"; tol.prior="dgamma(0.000001,0.000001)";
# #ar.prior="dnorm(0,1)"; tol.prior="dgamma(0.000001,0.000001)"; var.prior=NULL; sd.prior=NULL; mean.prior=ar.prior
# var.bugs<-function(y, ar.order=1, k=NULL, sim=FALSE, 
#                    mean.centre=FALSE, beg=ar.order+1,
#                    mean.prior=ar.prior, ar.prior="dnorm(0,1)", tolm.d.prior=0.001, tolm.nd.prior=0.001
#                    space=FALSE){
#   k=5
#   n<-dim(y)[1]
#   j<-dim(y)[2]
#   y<-matrix(y,n,j)
#   if(!is.null(k)){
#     y<-rbind(y,matrix(NA,k,j))
#   }
#   #k<-n-max(which(!is.na(y)))
#   if(beg<ar.order)
#     stop("The value of beg must be at least 1 greater than the number of lags")
# 
#   bug<-c("model{","")
#   #likelihood
#   lik<-c("#likelihood",
#          paste0("for(t in ",beg,":",n+k,"){"),
#          paste0("\ty[t] ~ dmnorm(y.mean[t,1:",j,"], isigma2[1:",j,",1:",j,"])"),
#          "}")
#   bug<-c(bug, lik)
#   #ymean
#   ymean<-c("#mean",
#            paste0("for(t in ",beg,":",n+k,"){"),
#            paste0("\ty.mean[t,",1:j,"] <- 0"))
#   if(ar.order==0 & mean.centre==TRUE) ymean[-(1:2)]<-paste0("\ty.mean[t,",1:j,"] <- phi0[",1:j,"]")
#   if(ar.order!=0 & mean.centre==FALSE){
#     temp<-paste0("\ty.mean[t,",1:j,"] <- ",paste0("phi",rep(1:ar.order,each=j),"[xx,",1:j,"]*y[t-",rep(1:ar.order,each=j),",",1:j,"]",collapse=" + "))
#     for(i in 1:j) temp[i] <- gsub("xx",i,temp[i])
#     ymean[-(1:2)]<-temp
#   }  
#   if(ar.order!=0 & mean.centre==TRUE){
#     temp<-paste0("\ty.mean[t,",1:j,"] <- phi0[",1:j,"] + ",paste0("phi",rep(1:ar.order,each=j),"[xx,",1:j,"]*(y[t-",rep(1:ar.order,each=j),",",1:j,"]-phi0[",1:j,"])",collapse=" + "))
#     for(i in 1:j) temp[i] <- gsub("xx",i,temp[i])
#     ymean[-(1:2)]<-temp
#   } 
#   bug<-c(bug, ymean, "}", "")
#   
# 
#   ####still to do....
#   ####not wishart prior could be diag(0.0001, 0.0001) as in camel example.
#   #prior
#   for(i in 1:j){
#     paste0("phi[",i,",",1:j,"] ~ ",ar.prior)
#   }
#   if(!is.null(tol.prior)){
#     ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
#                  paste0("isigma2 ~ ",tol.prior),
#                  "sigma <- pow(isigma2,-0.5)")
#   }
#   if(!is.null(var.prior)){
#     ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
#                  paste0("sigma2 ~ ",var.prior),
#                  "isigma2 <- pow(sigma2,-1)")
#   }
#   if(!is.null(sd.prior)){
#     ar.priors<-c(paste0("phi",0:ar.order," ~ ",ar.prior), 
#                  paste0("sigma ~ ",sd.prior),
#                  "isigma2 <- pow(sigma,-2)")
#   }
#   ar.priors<-ar.priors[-1]
#   if(mean.centre==TRUE)  ar.priors<-c(paste0("phi0 ~ ",mean.prior),ar.priors)
#   bug<-c(bug,"#priors",ar.priors,"")
#   
#   #forecast
#   forc<-NULL
#   if(k!=0){
#     forc<-c("#forecast",
#             paste("for(t in ",n+1,":",n+k,"){",sep=""),
#             "\ty.new[t] <- y[t]",
#             "}",
#             "")
#     bug<-c(bug,forc)
#   }  
#   
#   #simulation
#   ysim<-NULL
#   if(sim==TRUE){
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
#   #print.tsbugs(list(bug=bug))
#   
#   if(space==FALSE){
#     bug<-bug[-nchar(bug)!=0]
#     if(length(grep("#mean", bug))>0)
#       bug<-bug[-grep("#mean", bug)]
#   }
#   
#   p1<-grep("#likelihood",bug)
#   p2<-grep("#priors",bug)
#   if(k!=0 & sim==TRUE){
#     p3<-grep("#forecast",bug); p4<-grep("#simulation",bug)
#   }
#   if(k!=0 & sim==FALSE){
#     p3<-grep("#forecast",bug); p4<-length(bug)
#   }
#   if(k==0 & sim==TRUE){
#     p3<-grep("#simulation",bug); p4<-p3
#   } 
#   if(k==0 & sim==FALSE){
#     p3<-length(bug); p4<-p3
#   } 
#   p5<-length(bug)
#   
#   bug<-list(bug=bug,
#             data=list(y=y),
#             info=list(n=n,k=k,nh=n+k,
#                       args=mget(names(formals()),sys.frame(sys.nframe()))[-1],
#                       variance="MCV",
#                       likelihood=p1:(p2-1),
#                       priors=p2:(p3-1),
#                       forecast=NULL,
#                       simulation=NULL))
#   if(p3!=p4)  bug$info$forecast<-p3:(p4-1)
#   if(p4!=p5)  bug$info$simulation<-p4:(p5-1)
#   class(bug)<-"tsbugs"
#   return(bug)
# }
# 
# y[t,1] = c[1] + alpha[1,1] * y[t-1,1] + alpha[1,2] * y[t-1,2] + alpha[1,3] * y[t-1,3] + err1[t]
