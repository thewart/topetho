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
	int n;		//number of observations
	int K;		//number of clusters
	int mp;		//number of state variables
	int s[mp];	//number of states
	
	int Xp[n,sum(s)];	//state occupation counts in seconds
}

parameters {
	real<lower=0> lambda[K];	//unnormalized cluster weights
	real<lower=0> alpha_raw[K,sum(s)];	//unnormalized dirichlet parameters
	real<lower=0> A[K,mp];		//inverse sum of alphas
}

transformed parameters {
	real<lower=0> alpha[K,sum(s)];	//dirichlet parameters
	simplex[K] pi;			//cluster membership probabilities
	
	for (k in 1:K) {
		int pos;
		pos <- 1;
		pi[k] <- lambda[k] / sum(lambda);
		
		for (i in 1:mp) {
			real anorm;

			anorm <- sum(segment(alpha_raw[k],pos,s[i]));
			for (j in pos:(pos+s[i]-1))
				alpha[k,j] <- A[k,i] * (alpha_raw[k,j]/anorm);
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

			logp[k] <- log(pi[k]);
			for (j in 1:mp) {
				logp[k] <- logp[k] + dirmult_log(
					segment(Xp[i],pos,s[j]),s[j],
					segment(alpha[k],pos,s[j]),A[k,j]);
				pos <- pos + s[j];
			}
		}
		increment_log_prob(log_sum_exp(logp));
	}

	lambda ~ gamma(1.0/K,1.0);
	for (k in 1:K) {
		alpha_raw[k] ~ gamma(1,1);
		A[k] ~ cauchy(0,10);
	}
}
