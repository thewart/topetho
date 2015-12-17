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
init <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=40,chains=1,cores = 1)
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
ptetho <- c("Scratch","SelfGrm","Vigilnce","threat","Approach","Leave","noncontactAgg","InsptInf",
            "avoid","FearGrm","AffVoc","contactAgg","displace","Submit","GrmTerm","GrmPrsnt","ChinThrst")

Xc <- countprep(ptetho,bdat)
standat <- list(n=nrow(Xc),mc=length(ptetho),mp=length(statetho),
                s=ss,Xp=sdat,Xc=unlist(Xc),K=1,n0=1)
#with(standat,stan_rdump(ls(),file ="~/analysis/topetho/dump"))
#with(init,stan_rdump(ls(),file ="~/analysis/topetho/init"))

init <- stan(file="~/code/topetho/topetho.stan",data=standat,iter=100,chains=1,cores = 1,control=list(adapt_delta=0.6))
K <- 20
init <- list(theta_p_raw=summary(init)$summary[2:(sum(ss)+1),1] %>% 
               matrix(nrow=K,ncol=sum(ss),byrow = T),
             A=summary(init)$summary[(sum(ss)+2):(sum(ss)+1+length(statetho)),1] %>% 
               matrix(nrow=K,ncol=length(statetho),byrow=T),
             theta_c=summary(init)$summary[(sum(ss)+2+length(statetho)):(sum(ss)+1+length(statetho)+length(ptetho)),1] %>%
               matrix(nrow=K,ncol=length(ptetho),byrow=T),
             lambda=rep(1,K))
standat$K <- K
bigfit <- stan("~/code/topetho/topetho.stan",data = standat,chains = 2,cores = 2,
            init = list(init,init),iter = 200,control=list(adapt_delta=0.5,max_treedepth=8))
