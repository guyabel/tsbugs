##
##fn_tsbugs1.R: uv models
##fn_tsbugs2.R: general functions
##fn_tsbugs4.R: garch models
##fn_tsbugs4.R: mv models
##

print.tsbugs<-function(x, ...){
  cat(x$bug,sep="\n", ...)
}

inits<-function(bug, warn.mess=TRUE){
  if(class(bug)!="tsbugs")
    stop("bug must be object with tsbugs")
  in0<-nodes(bug, "prior")
  if(bug$info$variance=="RV") in0<-rbind(in0,nodes(bug, "likelihood"))
  in0<-in0[in0$stoc==1,]
  if(bug$info$variance!="RV")  in0<-in0[is.na(in0$beg) & is.na(in0$end),]
  if(bug$info$variance=="RV")  in0<-in0[in0$name!="y",]
  n<-dim(in0)[1]
  in1<-as.list(rep(0,n))
  names(in1)<-in0$name
  if(length(grep("sig",names(in1)))>0)  in1[[grep("sig",names(in1))]]<-0.5
  if(length(grep("tau",names(in1)))>0)  in1[[grep("tau",names(in1))]]<-0.5
  if(length(grep("psi",names(in1)))>0)  for(i in grep("psi",names(in1)))  in1[[i]]<-0.5
  if(length(grep("psi0.star",names(in1)))>0)  in1[[grep("psi0.star",names(in1))]]<-20
  if(bug$info$variance=="RV"){
    if(length(grep("lambda",names(in1)))>0)  
      in1[[grep("lambda",names(in1))]]<-1
    if(length(grep("epsilon",names(in1)))>0)  
      in1[[grep("epsilon",names(in1))]]<-0.05
    if(length(grep("delta",names(in1)))>0)  
      in1[[grep("delta",names(in1))]]<-c(rep(NA,in0[grep("delta",names(in1)),"beg"]),
                                         rep(0,in0[grep("delta",names(in1)),"end"]-in0[grep("delta",names(in1)),"beg"]))
    if(length(grep("beta",names(in1)))>0)
      in1[[grep("beta",names(in1))]]<-c(rep(NA,in0[grep("beta",names(in1)),"beg"]),
                                        rep(0,in0[grep("beta",names(in1)),"end"]-in0[grep("beta",names(in1)),"beg"]))
  }
  if(bug$info$variance=="GARCH"){
    if(length(grep("omega",names(in1)))>0)  for(i in grep("omega",names(in1)))  in1[[i]]<-0.5
    in1<-c(in1,list(h=c(rep(0.05,bug$info$args$vol.beg-1), rep(NA,bug$info$nh-bug$info$args$vol.beg+1))))
  }
  if(warn.mess==TRUE) print("guess attempt at initial values, might need to alter")
  return(in1)
}
# inits(garch0)


#bug=rv;part="likelihood"
nodes<-function(bug, part=NULL){
  if(class(bug)!="tsbugs")
    stop("bug must be object with tsbugs")
  if(!is.null(part))
    if(is.na(pmatch(part, names(bug$info))))
      stop("part must be a component of tsbug$info such as likelihood, priors, etc...")
  if(!is.null(part))
    bug<-bug$bug[bug$info[[pmatch(part, names(bug$info))]]]
  if(is.null(part))
    bug<-bug$bug
  
  nds<-gsub("\t","",bug)
  nds<-gsub("\\[t\\]","",nds)
  
  st.nds<-grep("~",nds)
  st<-nds[st.nds]
  st.dist<-gsub("(.*)~","",st)
  st.dist<-gsub("\\s","",st.dist)
  st<-gsub("~(.*)","",st)
  st<-gsub("\\s","",st)
  if(length(st)==0)  out<-NULL
  if(length(st)>0)  out<-data.frame(name=st,type="~",dt=st.dist,beg=NA,end=NA,stoc=1,id=st.nds)
  out$dist<-gsub("\\((.*)","",out$dt)
  out$param1<-gsub(",(.*)","",out$dt)
  out$param1<-gsub("(.*)\\(","",out$param1)
  out$param1[grep("[a-z]",out$param1)]<-NA
  #out$param2<-gsub("[a-z]","",out$param1)
  #out$param1<-gsub("^\\.","",out$param1)
  #out$param1<-gsub("\\.$","",out$param1)
  out$param1<-as.numeric(out$param1)
  out$param2<-gsub("(.*),","",out$dt)  
  out$param2<-gsub(")(.*)","",out$param2)
  out$param2[grep("[a-z]",out$param2)]<-NA
  #out$param2<-gsub("[a-z]","",out$param2)
  #out$param2<-gsub("^\\.","",out$param2)
  #out$param2<-gsub("\\.$","",out$param2)
  out$param2<-as.numeric(out$param2)
  
  det.nds<-grep("<-",nds)
  det<-nds[det.nds]
  det.trans<-gsub("(.*)<-","",det)
  det.trans<-gsub("\\s","",det.trans)
  det<-gsub("<-(.*)","",det)
  det<-gsub("\\s","",det)
  if(length(det)>0)  out<-rbind(out,data.frame(name=det,type="<-",dt=det.trans,beg=NA,end=NA,stoc=0,id=det.nds,
                                               dist=NA,param1=NA,param2=NA))
  
  st.loops<-grep("[0-9]:" ,nds)
  end.loops<-grep("}" ,nds)
  rg<-nds[st.loops]
  rg<-strsplit(gsub("[^[:digit:]:[:digit:]]","",rg),":")
  if(length(st.loops)>0){
    for(i in 1:dim(out)[1]){
      sel.node<-out$id[i]
      for(j in 1:length(st.loops)){
        if(sel.node>st.loops[j] & sel.node<end.loops[j]){
          out$beg[i]<-as.numeric(rg[[j]][1])
          out$end[i]<-as.numeric(rg[[j]][2])
        }
      }
    } 
  }
  out$name<-as.character(out$name)
  out
}
#nodes(rv, "likelihood")