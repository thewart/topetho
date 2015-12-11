functions {
	real ctmc_log(int[] ys, real[] yt, vector ps, 
		vector[] pt, real[] v, int len, real tlim) {

		real lp;
		lp <- 0;
		for (i in 1:len) {
			if (i==1) 			//initial state
				lp <- lp + log(ps[ys[i]]);
			else 				//transitioned-to state
				lp <- lp + log(pt[ys[i-1],ys[i]]) + 
					exponential_log(yt[i]-yt[i-1],v[ys[i-1]]);

			if (i==len)			//final state
				lp <- lp + exponential_ccdf_log(tlim-yt[i],v[ys[i]]);
			}
		return lp;
	}

	vector[] r2T(real[] r, int S) {
		vector[S] pt[S];
		for (s in 1:S) {
			for (s2 in 1:S) {
				if (s2==s) pt[s,s2] <- 0;
				else pt[s,s2] <- r[s2]/(sum(r)-r[s]);
				}
		}
		return pt;
	}
}
		

data {
	int n;				//number of observations
	int m;				//total number of state visitations
	int S;				//number of states
	int K;				//number of clusters
	real<lower=0> tlim;		//observation end time

	int<lower=1> len[n];		//per trial sequence length
	int<lower=0,upper=S> ys[m];	//state sequences
	real<lower=0,upper=tlim> yt[m];	//transition times (padded with initial zeros)

	real<lower=0> sigma0;		//prior scale for transition scales
	real<lower=0> n0;
}

parameters {
	real<lower=0> lambda[K];		//unnormalized cluster weights
	real<lower=0> theta[K,S];	//state transition scales
}

transformed parameters {
	simplex[K] pi;
	real<lower=0> r[K,S];		//state rates in
	real<lower=0> v[K,S];		//rates out
	simplex[S] ps[K];		//stationary distribution
	simplex[S] pt[K,S];		//transition matrix
	

	for (k in 1:K) {
		pi[k] <- lambda[k] / sum(lambda);
		for (s in 1:S)
			r[k,s] <- 1/theta[k,s];

		for (s in 1:S){
			ps[k,s] <- r[k,s]/sum(r[k]);
			v[k,s] <- sum(r[k]) - r[k,s];
		}	
		pt[k] <- r2T(r[k],S);
	}
}

model {
	int pos;
	pos <- 1;

	for (i in 1:n) {
		real logp[K];
		for (k in 1:K)
			logp[k] <- log(pi[k]) + ctmc_log( segment(ys,pos,len[i]),
				segment(yt,pos,len[i]),ps[k],pt[k],v[k],len[i],tlim);
		if (is_inf(log_sum_exp(logp)))
			reject("Observation ", i, " is impossible.");
		increment_log_prob( log_sum_exp(logp) );
		pos <- pos + len[i];
	}
	
	for (k in 1:K)
		theta[k] ~ cauchy(0,sigma0);
	lambda ~ gamma(n0/K,1);
}
