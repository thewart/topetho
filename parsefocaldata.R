derepeat <- function(behav) {
  n <- length(behav)
  x <- !logical(n)
  if (n>1) {
    for (i in 2:n) 
      if (behav[i-1]==behav[i]) x[i] <- F
  }
  return(x)
}

statvect <- function(time,behav,states) {
  out <- array(0,length(states))
  for (i in 1:length(behav))
    out[match(behav[i],states)] <- time[i]
  return(out)
}

statprep <- function(stat,dat,maxsec=630,nsec=10,ctmc=F) {
  sdat <- dat[Behavior %in% stat,.(Observation,Behavior,RelativeEventTime,Duration)]
  sdat <- sdat[,if ((length(Behavior)>1) & (RelativeEventTime[2] - RelativeEventTime[1]) < 2.6)
    .(Behavior[-1],RelativeEventTime[-1],Duration[-1])
    else .(Behavior,RelativeEventTime,Duration),by=Observation]
  
  if (!ctmc) {
    sdat[,V3:=pmin(V3,630-V2)]
    sdat <- sdat[,.(time=sum(round(V3/nsec))),by=c("Observation","V1")]
    sdat <- sdat[,statvect(time,V1,stat),by=Observation]
    out <- matrix(sdat$V1,ncol=length(stat),byrow=T)
  } else {
    sdat <- sdat[sdat[,derepeat(V1),by=Observation]$V1]
    sdat[,V2:=V2/(max(V2)+1)]
    setnames(sdat,c(2,3),c("Behavior","StartTime"))
    sdat[,V3:=NULL]
    out <- sdat
  }
  return(out)
}

countprep <- function(behav,dat) {
  
  pt <- dat[,length(EventName),by=c("Observation","Behavior")]
  pt <- dcast.data.table(pt,Observation ~ Behavior,fill=0)
  pt <- pt[,which(names(pt) %in% behav),with=F]
  return(as.matrix(pt))
}