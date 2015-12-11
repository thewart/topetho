derepeat <- function(behav) {
  n <- length(behav)
  x <- !logical(n)
  if (n>1) {
    for (i in 2:n) 
      if (behav[i-1]==behav[i]) x[i] <- F
  }
  return(x)
}

actvect <- function(time,behav,states) {
  out <- array(0,length(states))
  for (i in 1:length(behav))
    out[match(behav[i],states)] <- time[i]
  return(out)
}

library(data.table)
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv"))
                           
bdat <- fread(fpath)
bdat <- copy(bdat[!overtime & !BadObs])

actetho <- c("Rest","GroomGIVE","GroomGET","Feed","Travel")
activ <- bdat[Behavior %in% actetho,.(Observation,Behavior,RelativeEventTime,Duration)]
activ <- activ[,if ((length(Behavior)>1) & (RelativeEventTime[2] - RelativeEventTime[1]) < 2.6)
  .(Behavior[-1],RelativeEventTime[-1],Duration[-1])
  else .(Behavior,RelativeEventTime,Duration),by=Observation]
activ <- activ[activ[,derepeat(V1),by=Observation]$V1]
activ[,V3:=pmin(V3,630-V2)]
#activ$V2 <- activ[,V2/(max(V2)+1)]
#standat <- list(n=activ[,length(unique(Observation))], m=nrow(activ),
#                S=length(actetho),tlim=1,len=activ[,length(V1),by=Observation]$V1,
#                ys=as.numeric(factor(activ$V1)),yt=activ$V2,sigma0=2.5,K=10,n0=1)
# 
# init <- with(standat,list(
#   lambda=rep(n0/K,K),theta=matrix(unlist(activ[,mean(V3)/600,by=V1][,2,with=F]),K,S,byrow = T)))
# ctmcfit <- stan(file = "~/code/topetho/ctmc_clust_full.stan",data = standat,pars = c("pi","Q"),
#             chains = 2,iter = 300,warmup = 200,control = list(adapt_delta=0.65))

activ <- activ[,.(time=sum(round(V3/10))),by=c("Observation","V1")]
activ <- activ[,actvect(time,V1,actetho),by=Observation]
activ <- matrix(activ$V1,ncol=length(actetho),byrow=T)

standat <- list(s=length(actetho),n=nrow(activ),K=1,Xp=activ)
init <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=100,
             chains=1,cores = 1)
init <- list(alpha_raw=matrix(summary(init)$summary[2:6,1],nrow=20,ncol=5,byrow = T),
             Ainv=rep(summary(init)$summary[7,1],20),
             lambda=rep(1,20))
standat <- list(s=length(actetho),n=nrow(activ),K=20,Xp=activ)
dirfit_decs <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=300,
            init = list(init,init),chains=2,control=list(adapt_delta=0.6),cores = 2)

library(rstan)
rstan_options(auto_write = TRUE)

ptetho <- c("Vigilnce","threat","Approach","Leave","noncontactAgg",
            "avoid","FearGrm","AffVoc","contactAgg","displace","Submit")
pt <- bdat[,length(EventName),by=c("Observation","Behavior")]
pt <- dcast.data.table(pt,Observation ~ Behavior,fill=0)
pt <- pt[,which(names(pt) %in% ptetho),with=F]

Xp <- as.matrix(activ[,-1,with=F])
Xc <- as.matrix(pt)
#Xc <- apply(Xc,2,function(x) x/sd(x))

standat <- list(n=nrow(Xc),mc=length(ptetho),
                #mp=length(actetho),Xp=unlist(Xp),
                Xc=unlist(Xc),K=10,n0=1)
with(standat,stan_rdump(ls(),file ="~/analysis/topetho/dump"))

init <- with(standat,list(lambda=rep(n0/K,K),theta_c=matrix(colMeans(Xc),K,mc,byrow = T)))
with(init,stan_rdump(ls(),file ="~/analysis/topetho/init"))

foot <- stan("~/code/topetho/topetho_alt.stan",data = standat,chains = 1,
            init = list(init),iter = 10,warmup = 10)
