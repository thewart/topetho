functions {
	real ctmc_log(int[] ys, real[] yt, vector ps, vector[] pt, real[] v, int len, real tlim) {
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
	real<lower=0> tlim;		//observation end time

	int<lower=1> len[n];		//per trial sequence length
	int<lower=0,upper=S> ys[m];	//state sequences
	real<lower=0,upper=tlim> yt[m];	//transition times (padded with initial zeros)

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
	}	
	pt <- r2T(r,S);
}

model {
	int pos;
	pos <- 1;

	for (i in 1:n) {
		increment_log_prob( ctmc_log( segment(ys,pos,len[i]),
				segment(yt,pos,len[i]),ps,pt,v,len[i],tlim) );
		pos <- pos + len[i];
	}

	theta ~ cauchy(0,sigma0);
}
