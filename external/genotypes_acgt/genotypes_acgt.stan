//
// Model of strain genotypes at varying abundances
//

data {
	int <lower=0> U;         // number of SNV sites
	int <lower=0> T;         // number of time points
	int <lower=0> S;         // max number of strains
	int <lower=0> nota[U*T]; // vector of counts of A  
	int <lower=0> notc[U*T]; // vector of counts of C
	int <lower=0> notg[U*T]; // vector of counts of G  
	int <lower=0> nott[U*T]; // vector of counts of T  
}

transformed data {
	int a[U,T]; // count of A matrix, site x timepoint
	int c[U,T]; // count of C matrix, site x timepoint
	int g[U,T]; // count of G matrix, site x timepoint
	int t[U,T]; // count of T matrix, site x timepoint
	for(x in 1:U*T){
		a[((x-1)/T)+1,((x-1)%T)+1] <- nota[x];
		c[((x-1)/T)+1,((x-1)%T)+1] <- notc[x];
		g[((x-1)/T)+1,((x-1)%T)+1] <- notg[x];
		t[((x-1)/T)+1,((x-1)%T)+1] <- nott[x];
	}
}

parameters {	
	vector <lower=0>[S] depth[T]; // strain depth of coverage at each time point
    simplex[4] strain[U,S];
}

model {
	// set the priors
	for( i in 1:T ) {
		depth[i] ~ uniform(0,200); // uniformly distributed strain abundances
	}
	// calculate log likelihood
	for( u in 1:U) {
		for( i in 1:T ) {
            vector[S] strain_a;
            vector[S] strain_c;
            vector[S] strain_g;
            vector[S] strain_t;
            for( s in 1:S ){
                strain_a[s] <- strain[u,s,1];
                strain_c[s] <- strain[u,s,2];
                strain_g[s] <- strain[u,s,3];
                strain_t[s] <- strain[u,s,4];
            }
			// poisson counts of A, C, G, T
			increment_log_prob( poisson_log( a[u,i], depth[i]'*strain_a) );
			increment_log_prob( poisson_log( c[u,i], depth[i]'*strain_c) );
			increment_log_prob( poisson_log( g[u,i], depth[i]'*strain_g) );
			increment_log_prob( poisson_log( t[u,i], depth[i]'*strain_t) );
		}
	}
}

