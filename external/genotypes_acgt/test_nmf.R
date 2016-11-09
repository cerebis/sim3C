#!/usr/bin/env Rscript
library(rstan)
#
# define the number of strains, samples, and sites
strains <- 6
samples <- 8
sites <- 50

maxcov <- 100
#
# choose random values for strain abundances at all timepoints
test_abunds<-matrix(runif(strains*samples, 0, maxcov), strains, samples)

#
# choose random values for strain genotypes, ensuring at least 1 mutation at all positions
primary <- round(runif(sites, 0.501, strains+0.49))
test_sites <- matrix(rep(0,strains*sites),strains,sites)
for(i in 1:sites){
    test_sites[primary[i],i]<-1
}
zlen <- sum(test_sites==0)
test_sites[test_sites==0] <- round(runif(zlen,0,1))

#
# generate the matrix of mutation frequency (ground truth)
mixed<-round(t(test_abunds)%*%test_sites)
total<-round(t(test_abunds)%*%matrix(rep(1,strains*sites),strains,sites))
#
# generate a Poisson distributed number of site observations with rpois()
observation_counts<-matrix(rpois(samples*sites, total),samples,sites)
#
# A binomially distributed fraction of the observations are the mutant
observed_mut<-matrix(rbinom(samples*sites,observation_counts,mixed/total),samples,sites)

bestorder <- function(x,y) { 
    library(combinat)
    ppp<-permn(seq(1,nrow(x)))
    bestgeno <- ""
    bestscore <- 99999999
    for(i in 1:length(ppp)){
        cs <- x[ppp[[i]],]-y
        cs <- cs * cs
        cs <- sqrt(sum(cs)/length(cs))
        if(cs < bestscore){
            bestgeno <- x[ppp[[i]],]
            bestscore <- cs
        }
    }
    print(paste("min RMSE ",bestscore))
    bestgeno
}


#
# define the stan model as a text string
stancode <- "
data {
	int <lower=0> U;         // number of SNV sites
	int <lower=0> T;         // number of time points
	int <lower=0> S;         // max number of strains
	int <lower=0> observations[U,T]; // matrix of read counts at each site & timepoint
	int <lower=0> mutations[U,T]; // matrix of observed mutation counts at each site & timepoint
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
	for( i in 1:U ) {
    	for( j in 1:S ) {
		    geno[i][j] ~ beta(0.1,0.1); // prior concentrated on 0 & 1 (poor man's approximation to what should be binary variable)
        }
	}
	// calculate log likelihood
	for( u in 1:U) {
		for( i in 1:T ) {
            vector[S] dnorm;
			// Poisson distributed number of reads
			increment_log_prob( poisson_log( observations[u,i], sum(depth[i])) );
            // Binomial distributed number of SNV observations
            dnorm <- depth[i] / sum(depth[i]);
			increment_log_prob( binomial_log( mutations[u,i], observations[u,i], dnorm'*geno[u]) );
		}
	}
}
"


sm <- stan_model(model_name="genotypes",model_code=stancode)
dat <- list(U = sites, T = samples, S = strains, observations = t(observation_counts), mutations=t(observed_mut));
vbfit <- vb(sm, data = dat, tol_rel_obj=0.001, elbo_samples=200, iter=20000, sample_file = 'geno.csv', algorithm="meanfield")


# print(vbfit)

"True depths"
test_abunds
""

ddest<-matrix(get_posterior_mean(vbfit,pars="depth"),strains,samples)
"Posterior mean estimates"
ddest

# now compute a factorization with the NMF package
library("NMF")
nnn<-nmf(observed_mut / observation_counts,strains,nrun=500)

ccc<-coef(nnn)
ccc<-bestorder(ccc,test_sites)

geno<-matrix(get_posterior_mean(vbfit,pars="geno"),strains,sites)
geno <- bestorder(geno, test_sites)

x11(width=8)
par(mfrow=c(3,1),mar=c(2,1,2,1))
image(z=t(test_sites),,col=gray((32:0)/32),main="Truth",xaxt="n")
image(z=t(1-geno),col=gray((0:32)/32),main="MFVBBPNMF",xaxt="n")
image(z=t(ccc),,col=gray((32:0)/32),main="R NMF",xaxt="n")




dat <- list(U = sites, T = samples, S = strains-2, observations = t(observation_counts), mutations=t(observed_mut));
vbfit <- vb(sm, data = dat, tol_rel_obj=0.004, elbo_samples=200, iter=20000, sample_file = 'geno.csv', algorithm="meanfield")

dat <- list(U = sites, T = samples, S = strains-1, observations = t(observation_counts), mutations=t(observed_mut));
vbfit <- vb(sm, data = dat, tol_rel_obj=0.004, elbo_samples=200, iter=20000, sample_file = 'geno.csv', algorithm="meanfield")

dat <- list(U = sites, T = samples, S = strains, observations = t(observation_counts), mutations=t(observed_mut));
vbfit <- vb(sm, data = dat, tol_rel_obj=0.004, elbo_samples=200, iter=20000, sample_file = 'geno.csv', algorithm="meanfield")


dat <- list(U = sites, T = samples, S = strains+1, observations = t(observation_counts), mutations=t(observed_mut));
vbfit <- vb(sm, data = dat, tol_rel_obj=0.004, elbo_samples=200, iter=20000, sample_file = 'geno.csv', algorithm="meanfield")

