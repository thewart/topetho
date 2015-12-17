clusts <- c(2,3,5,6,7)
moo <- list()
nsamp <- with(attributes(fit)$stan_args[[1]],iter - warmup)
for (i in 1:length(clusts))
{
  moo[[i]] <- extract(fit,pars=paste0("r[",clusts[i],",",1:6,"]"))
  names(moo[[i]]) <- actetho
  moo[[i]] <- melt(as.data.table(moo[[i]]))
  moo[[i]]$cluster <- as.character(clusts[i])
  moo[[i]]$iter <- 1:nsamp
}
moo <- as.data.table(plyr::ldply(moo))
moo <- moo[!(variable %in% c("OutCorral","nopascon","NoGrmInf"))]
moo$cluster <-ordered(moo$cluster,levels=as.character(clusts))
ggplot(moo,aes(x=value,colour=cluster)) + geom_density() + facet_wrap(~variable,scales = "free")

foo <- moo[,.(cluster,rank=as.vector(t(apply(matrix(value,ncol = length(unique(cluster)))
                  ,1,function(x) rank(x)/length(x))))),by=variable]
foo <- foo[,.(rank=median(rank),lb=quantile(rank,0.05),ub=quantile(rank,0.95)),by=c("variable","cluster")]

ggplot(foo,aes(y=variable,x=rank)) + geom_point() +
  geom_errorbarh(aes(xmax=ub,xmin=lb),height=0)+ facet_grid(~cluster)


#comparing chains
fit <- dirfit_idecs
piplt <- rbind(melt(as.data.table(extract(fit,pars="pi",permuted=F)[,1,]))[,chain:=1],
             melt(as.data.table(extract(fit,pars="pi",permuted=F)[,2,]))[,chain:=2])
ggplot(piplt[value>0.04],aes(x=value,color=variable)) + geom_density() + 
  facet_grid(chain~.,scales = "free_y") + theme_classic() + theme(legend.position="none")

#dir mult alphas
moo <- clustercontent(bigfit,"theta_p",statetho,"theta_c",colnames(Xc),exclude = c("OutCorral","nopascon","NoGrmInf"),chain=1)
ggplot(moo[cluster %in% moo[,mean(weight)>0.1,by=cluster][V1==T,cluster]],
            aes(x=value,colour=cluster)) + geom_density() + facet_wrap(~variable,scales = "free")

#pull one chain
moo <- extract(ctmcfit,pars="pi",inc_warmup=T,permute=F)[,2,]
moo <- melt(as.data.table(moo))[,iter:=1:300]
ggplot(moo,aes(x=iter,y=value)) + geom_line() + facet_wrap(~variable,scale="free")
