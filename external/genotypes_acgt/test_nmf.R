#!/usr/bin/env Rscript
library(rstan)
#
# define the number of strains, samples, and sites
strains <- 2
samples <- 12
sites <- 300

#
# choose random values for strain abundances at all timepoints
test_abunds<-matrix(runif(strains*samples, 0, 100), strains, samples)

#
# choose random values for strain genotypes, ensuring at least 1 mutation at all positions
primary <- round(runif(sites, 1, strains))
test_sites <- matrix(rep(0,strains*sites),strains,sites)
for(i in 1:sites){
    test_sites[primary[i],i]<-1
}
zlen <- sum(test_sites==0)
test_sites[test_sites==0] <- round(runif(zlen,0,1))

#
# generate the matrix of mutation frequency
mixed<-round(t(test_abunds)%*%test_sites)

#
# generate observations of mutation frequency that are Poisson distributed with rpois()
observed<-matrix(rpois(samples*sites,mixed),samples,sites)

#
# define the stan model as a text string
stancode <- "
data {
	int <lower=0> U;         // number of SNV sites
	int <lower=0> T;         // number of time points
	int <lower=0> S;         // max number of strains
	int <lower=0> mixed[U,T]; // matrix of observed mutation counts at each site & timepoint
}

parameters {	
	vector <lower=0>[S] depth[T]; // strain depth of coverage at each time point
    vector<lower=0,upper=1>[S] geno[U];
}

model {
	// set the priors
	for( i in 1:T ) {
		depth[i] ~ uniform(0,500); // uniformly distributed strain abundances
	}
	// calculate log likelihood
	for( u in 1:U) {
		for( i in 1:T ) {
			// poisson prob
			increment_log_prob( poisson_log( mixed[u,i], depth[i]'*geno[u]) );
		}
	}
}
"


dat <- list(U = sites, T = samples, S = strains, mixed = t(observed));
sm <- stan_model(model_name="genotypes",model_code=stancode)
vbfit <- vb(sm, data = dat, sample_file = 'geno.csv', algorithm="fullrank")

"True depths"
test_abunds
""

print(vbfit)


