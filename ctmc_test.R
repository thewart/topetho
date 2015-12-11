#markov process simulation
library(markovchain)
library(rstan)
rstan_options(auto_write = TRUE)

# rates <- c(1,3,5)
# Qd <- c(sum(-rates[2:3]),sum(-rates[c(1,3)]),sum(-rates[1:2]))
# Q <- matrix(c(Qd[1],rates[2:3],rates[1],Qd[2],rates[3],rates[1:2],Qd[3]),
            # nrow=3,ncol=3,byrow=T)
Q <- matrix(c(-5,3,2,2,-6,4,.5,2,-2.5),3,byrow = T)
guh <- new("ctmc",byrow = T, generator = Q)
n <- 1000
tmax <- 1
tlst <- slst <- list()
for (i in 1:n)
{
  juh <- rctmc(Inf,guh,initDist = steadyStates(guh),T = tmax)
  tlst[[i]] <- juh[[2]]
  slst[[i]] <- as.numeric(juh[[1]])
}


standat <- list(n=n,m=sum(sapply(slst,length)),S=length(rates),tmax=tmax,
                tnum=sapply(slst,length),stseq=unlist(slst),tnt=unlist(tlst),sigma0=1)

moo <- stan(file = "~/code/topetho/ctmc_test.stan",data = standat,
            chains = 1,iter = 200,warmup = 100)

standat <- list(n=n,m=sum(sapply(slst,length)),S=length(rates),tlim=tmax,
                len=sapply(slst,length),ys=unlist(slst),yt=unlist(tlst),sigma0=1)
moo <- stan(file = "~/code/topetho/ctmc_full.stan",data = standat,
            chains = 1,iter = 200,warmup = 100)

standat$K <- 5
standat$n0 <- 1
moo <- stan(file = "~/code/topetho/ctmc_clust_full.stan",data = standat,
            chains = 1,iter = 200,warmup = 100,control=list(adapt_delta=0.65))
