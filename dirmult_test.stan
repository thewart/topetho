functions {
	real dirmult_log(int[] y, int s, real[] alpha, real A) {
		int N;
		real lp;
		N <- sum(y);
		lp <- log(N) + lbeta(A,N);
		for (i in 1:s) 
			if (y[i] > 0)
				lp <- lp - log(y[i]) - lbeta(alpha[i],y[i]);
		return lp;
	}
}
data {
	int s;		//number of states
	int n;		//number of observations
	int K;		//number of clusters
	
	int Xp[n,s];	//state occupation counts in seconds
}

parameters {
	real<lower=0> lambda[K];		//unnormalized cluster weights
	real<lower=0> alpha_raw[K,s];	//unnormalized dirichlet parameters
	real<lower=0> A[K];		//inverse sum of alphas
}

transformed parameters {
	real<lower=0> alpha[K,s];		//dirichlet parameters
	simplex[K] pi;		//cluster membership probabilities
	
	for (k in 1:K) {
		pi[k] <- lambda[k] / sum(lambda);
		for (i in 1:s) 
			alpha[k,i] <- A[k] * (alpha_raw[k,i]/sum(alpha_raw[k]));
	}
}

model {

	for (i in 1:n) {
		real logp[K];
		for (k in 1:K) 
			logp[k] <- log(pi[k]) + dirmult_log(Xp[i],s,alpha[k],A[k]);
		increment_log_prob(log_sum_exp(logp));
	}

	lambda ~ gamma(1.0/K,1.0);
	for (k in 1:K)
		alpha_raw[k] ~ gamma(1,1);
	A ~ cauchy(0,10);
}
