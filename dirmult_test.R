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
