functions {
	real ctmc_log(int[] ys, real[] yt, vector ps, vector[] pt, vector v, int len, real tlim) {
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

	matrix t2Q(real[,] theta, int S) {
		matrix[S,S] Q;
		for (s1 in 1:S) {
			for (s2 in 1:S) {
				if (s1>s2) Q[s1,s2] <- theta[s1,s2];
				else if (s1<s2) Q[s1,s2] <- theta[s1,s2-1];
				else Q[s1,s2] <- 0;
			}
		}
		for (s in 1:S) Q[s,s] <- -sum(Q[s]);
		return Q;
	}
	
	vector Q2p(matrix Q, int S) {
		matrix[S,S] A;
		vector[S] b;
		A <- append_row( block(Q,1,1,S,S-1)', rep_row_vector(1,S) );
		b <- rep_vector(0,S);
		b[S] <- 1;
		return A \ b;
	}

	vector[] Q2T(matrix Q, int S) {
		vector[S] pt[S];
		for (s1 in 1:S) {
			for (s2 in 1:S) {
				if (s2==s1) pt[s1,s2] <- 0;
				else pt[s1,s2] <- -Q[s1,s2]/Q[s1,s1];
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
	real<lower=0> theta[S,S-1];	//state transition scales
}

transformed parameters {
	matrix[S,S] Q;			//generator matrix
	simplex[S] ps;			//stationary distribution
	simplex[S] pt[S];		//transition matrix

	Q <- t2Q(theta,S);
	ps <- Q2p(Q,S);
	pt <- Q2T(Q,S);
}

model {
	int pos;
	pos <- 1;

	for (i in 1:n) {
		increment_log_prob( ctmc_log( segment(ys,pos,len[i]),
				segment(yt,pos,len[i]),ps,pt,-diagonal(Q),len[i],tlim) );
		pos <- pos + len[i];
	}

	for (s in 1:S)	theta[s] ~ cauchy(0,sigma0);
}
