#load scripts and data
library(data.table)
library(rstan)
rstan_options(auto_write = TRUE)
source("~/code/topetho/parsefocaldata.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv"))
                           
bdat <- fread(fpath)
bdat <- copy(bdat[!overtime & !BadObs])

#construct state data
statetho <- list(act=c("Rest","GroomGIVE","GroomGET","Feed","Travel"),
                 corr=c("OutCorral","InCorral"),
                 passetho=c("nopascon","passcont"),
                 inf=c("NoGrmInf","GromInf"))
ss <- sapply(statetho,length)
sdat <- lapply(statetho,statprep,dat=bdat)
sdat <- unlist(lapply(sdat,as.vector)) %>% matrix(ncol=sum(sapply(statetho,length)))

#fit model
standat <- list(mp=length(statetho),s=ss,n=nrow(sdat),K=1,Xp=sdat)
init <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=40,
             chains=1,cores = 1)
K <- 15;
init <- list(alpha_raw=summary(init)$summary[2:(sum(ss)+1),1] %>% 
               matrix(nrow=K,ncol=sum(ss),byrow = T),
             A=summary(init)$summary[(sum(ss)+2):(sum(ss)+1+length(statetho)),1] %>% 
               matrix(nrow=K,ncol=length(statetho),byrow=T),
             lambda=rep(1,K))
standat$K <- K;
dirfit <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=200,
            init = list(init,init),chains=2,control=list(adapt_delta=0.6),cores = 2)

#construct point data
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
