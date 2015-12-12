library(gtools)
n <- 1000
a <- c(2,5,1,0.5)
p <- rdirichlet(n,a)
X <- matrix(nrow=n,ncol=length(a))
for (i in 1:n)
  X[i,] <- rmultinom(1,20,p[i,])

library(rstan)
rstan_options(auto_write = TRUE)
standat <- list(s=length(a),n=n,Xp=X,K=10)
moo <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=200,warmup=100,
            chains=1,control=list(adapt_delta=0.65))

a2 <- c(1,0.25,0.75)
p <- rdirichlet(n,a2)
X2 <- matrix(nrow=n,ncol=length(a2))
for (i in 1:n)
  X2[i,] <- rmultinom(1,10,p[i,])
X <- cbind(X,X2)

standat <- list(mp=2,s=c(length(a),length(a2)),n=n,K=1,Xp=X)
moo <- stan(file="~/code/topetho/dirmult_test.stan",data=standat,iter=200,warmup=100,
            chains=1,control=list(adapt_delta=0.65))
