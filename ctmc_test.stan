data {
	int n;				//number of observations
	int m;				//total number of state visitations
	int S;				//number of states
	real<lower=0> tmax;		//observation end time

	int<lower=1> tnum[n];		//per trial sequence length
	int<lower=0,upper=S> stseq[m];	//state sequences
	real<lower=0,upper=1> tnt[m]; 	//transition times (padded with initial zeros)

	real<lower=0> sigma0;		//prior scale for transition scales
}

parameters {
	real<lower=0> theta[S];		//state transition scales
}

transformed parameters {
	real<lower=0> r[S];		//state rates in
	real<lower=0> v[S];		//rates out
	simplex[S] ps;			//stationary distribution
	simplex[S] pt[S];		//transition matrix

	for (s in 1:S)
		r[s] <- 1/theta[s];

	for (s in 1:S){
		ps[s] <- r[s]/sum(r);
		v[s] <- sum(r) - r[s];
		
		for (s2 in 1:S) {
			if (s2==s) pt[s,s2] <- 0;
			else pt[s,s2] <- r[s2]/(sum(r)-r[s]);
		}
	}

}

model {
	int pos;
	pos <- 1;

	for (i in 1:n) {
		for (j in pos:(pos+tnum[i]-1)) {
			if (j==pos) 			//initial state
				increment_log_prob( log(ps[stseq[j]]) );
			else {
				increment_log_prob( log(pt[stseq[j-1],stseq[j]]) );
				tnt[j]-tnt[j-1] ~ exponential(v[stseq[j-1]]);
			}
			if (j==pos+tnum[i]-1)	//final state
				increment_log_prob( 
					exponential_ccdf_log(tmax-tnt[j],v[stseq[j]]) );
		}
		
		pos <- pos + tnum[i];
	}

	theta ~ normal(0,sigma0);
}
