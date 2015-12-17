clustercontent <- function(fit,stvar="theta_p",stname,ptvar="theta_c",ptname,weight="pi",incthresh=0.001,chain="all",exclude=c(),normalpha=T) {
  stln <- c(0,cumsum(sapply(stname,length)))
  if (chain=="all")
    chain <- 1:ncol(extract(fit,"lp__",permute=F))
  
  out <- list()
  l <- 0
  for (k in chain) {
    pi <- extract(fit,pars=weight,permuted=F)[,k,]
    clusts <- which(colMeans(pi)>incthresh)
    for (i in 1:length(clusts)) {
      l <- l+1
      out[[l]] <- extract(fit,pars=paste0(stvar,"[",clusts[i],",",1:max(stln),"]"),permute=F)[,k,]
      if (normalpha) {
        for (j in 2:length(stln)) {
          ind <- (stln[j-1]+1):stln[j]
          out[[l]][,ind] <- out[[l]][,ind]/rowSums(out[[l]][,ind])
        }
      }
      out[[l]] <- cbind(out[[l]],extract(fit,pars=paste0(ptvar,"[",clusts[i],",",1:length(ptname),"]"),permute=F)[,k,])
      colnames(out[[l]]) <- c(unlist(stname),ptname)
      out[[l]] <- melt(as.data.table(out[[l]]))
      out[[l]]$chain <- k
      out[[l]]$weight <- pi[,clusts[i]]
      out[[l]]$cluster <- paste(clusts[i],"(",format(mean(out[[l]]$weight),digits=2),")")
      out[[l]]$iter <- 1:nrow(pi)
    }
  }
  out <- as.data.table(plyr::ldply(out))
  out <- out[!(variable %in% exclude)]
  return(out)
}