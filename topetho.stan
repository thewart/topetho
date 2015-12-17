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
	int n; 					//number of observations
	int mp; 				//number of proportional behaviors
	int mc;					//number of count behaviors
	int s[mp];				//number of states prop behavior

	int<lower=0> Xp[n,sum(s)];		//data for proportional behaviors
	int<lower=0> Xc[n,mc];			//data of count behaviors

	int K;					//maximum number of clusters
	real n0;				//prior sample size
}

parameters {
	real<lower=0> lambda[K];		//unnormalized cluster weights

	real<lower=0> theta_p_raw[K,sum(s)];	//unnormalized dirichlet parameters
	real<lower=0> A[K,mp];			//sum of dir params
	
	real<lower=0> theta_c[K,mc];		//means for count data	
}

transformed parameters {
	simplex[K] pi;				//cluster membership probabilities
	real<lower=0> theta_p[K,sum(s)];	//dirichlet parameters
	
	for (k in 1:K) {
		int pos;
		pos <- 1;
		pi[k] <- lambda[k] / sum(lambda);
		
		for (i in 1:mp) {
			real anorm;

			anorm <- sum(segment(theta_p_raw[k],pos,s[i]));
			for (j in pos:(pos+s[i]-1))
				theta_p[k,j] <- A[k,i] * (theta_p_raw[k,j]/anorm);
			pos <- pos + s[i];
		}
	}
}

model {
	for (i in 1:n) {
		real logp[K];
		for (k in 1:K) {
			int pos;
			pos <- 1;

			logp[k] <- log(pi[k]) + poisson_log(Xc[i],theta_c[k]);
			for (j in 1:mp) {
				logp[k] <- logp[k] + dirmult_log(
					segment(Xp[i],pos,s[j]),s[j],
					segment(theta_p[k],pos,s[j]),A[k,j]);
				pos <- pos + s[j];
			}
		}
		increment_log_prob(log_sum_exp(logp));
	}
	
	lambda ~ gamma(n0/K,1);
	for (k in 1:K) {
		theta_p_raw[k] ~ gamma(1,1);
		theta_c[k] ~ cauchy(0,10);
		A[k] ~ cauchy(0,10);
	}
}
