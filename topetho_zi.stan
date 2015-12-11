data {
	int n; 		//number of observations
	int mp; 	//number of proportional behaviors
	int mb;		//number of binary behaviors
	int mc;		//number of count behaviors

	simplex[mp] Xp[n];		//data for proportional behaviors
	int<lower=0,upper=1> Xb[n,mb];	//data of binary behaviors
	int<lower=0> Xc[n,mc];		//data of count behaviors

	int K;		//maximum number of clusters
	real n0;	//prior sample size
}

parameters {
	real<lower=0> lambda[K];		//unnormalized cluster weights

	real<lower=1,upper=1> p0c[K,mc];	//zero weight for counts
	vector<lower=0>[mp] theta_p[K];		//dirichlet params for prop behavs
	real<lower=0,upper=1> theta_b[K,mb]	//probabilities for binary behavs
	real<lower=0> theta_c[K,mc]		//means for count behavs
}

transformed parameters {
	simplex[K] pi;		//cluster membership probabilities
	
	for (k in 1:K)
		pi[k] <- lambda[k] / sum(lambda);
}

model {
	
	for (i in 1:n)
	{
		real logp[K];			//log probs for each behavior type
		real logpc[K];

		for (k in 1:K)
		{
			logpc[k] <- 0;
			for (j in 1:mc)
				logpc[k] <- logpc[k] + logsumexp(
					log1m(p0c[k,j]) + poisson_log(Xc[i,j],theta_c[k]),
					log(p0c[k,j]) + ifelse(Xc[i,j]==0,log(1),log(0)));

			logp[k] <- log(pi[k]) + dirichlet_log(Xp[i],theta_p[k]) + 
				binomial_log(Xb[i],theta_b[k]) + logpc[k];
		}
		
		increment_log_prob(log_sum_exp(logp));
	}
	
	lambda ~ gamma(n0/K,1);
	//theta_b ~ uniform(0,1);
	for (k in 1:K)
	{
		theta_p[k] ~ cauchy(0,1);
		theta_c[k] ~ cauchy(0,5);
	}	
}
