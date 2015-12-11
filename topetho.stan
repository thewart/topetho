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
	int n; 		//number of observations
	int mp; 	//number of proportional behaviors
	int mc;		//number of count behaviors

	int<lower=0> Xp[n,mp];		//data for proportional behaviors
	int<lower=0> Xc[n,mc];		//data of count behaviors

	int K;		//maximum number of clusters
	real n0;	//prior sample size
}

parameters {
	real<lower=0> lambda[K];		//unnormalized cluster weights

	real<lower=0> theta_p_raw[K,s];		//unnormalized dirichlet parameters
	real<lower=0> theta_p_scale[K];		//inverse sum of alphas
	
	real<lower=0> theta_c[K,mc];			//means for count data	
}

transformed parameters {
	simplex[K] pi;			//cluster membership probabilities
	real<lower=0> theta_p[K,s];	//dirichlet parameters
	
	for (k in 1:K)
		pi[k] <- lambda[k] / sum(lambda);
		for (i in 1:s) 
			theta_p[k,i] <- (1/theta_p_scale[k]) * 
				(theta_p_raw[k,i]/sum(theta_p_raw[k]));
	}
}

model {
	
	for (i in 1:n) {
		real logp[K];
		for (k in 1:K)
			logp[k] <- log(pi[k]) + poisson_log(Xc[i],theta_c[k]) +
			dirmult_log(Xp[i],mp,theta_p[k],1/theta_p_scale[k]);
		
		increment_log_prob(log_sum_exp(logp));
	}
	
	lambda ~ gamma(n0/K,1);
	for (k in 1:K) {
		theta_p_raw[k] ~ gamma(1,1);
		theta_c[k] ~ cauchy(0,10);
	}
	theta_p_scale ~ cauchy(0,10);
}
