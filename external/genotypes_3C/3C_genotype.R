#!/usr/bin/env Rscript
# (c) 2016 Aaron Darling
# 
library(rstan)


# formulation of genotype resolution from noisy 3C measurements of SNVs as
# a graph factorization problem
genotype3C <- "
data {
    int V; # number of sites (nodes in graph)
    int L; # number of links among sites (edges in graph)
    int linkcounts[L,4]; # counts of SNV co-observation. Genotypes are binary variables: mutant or not. Four possible combinations among two sites: 00,01,10,11 
    int linksites[L,2]; # the indices of the two sites involved in each link
    int K; # number of components to factorize
}

parameters {
    real<lower=0,upper=1> genotype[K,V]; # the genotype of the K strains. Ideally this would be a binary (or base 4) variable
    simplex[K] abundance;  # the abundance of each strain
}

model {
    for(k in 1:K){
        for(v in 1:V){
            genotype[k,v] ~ beta(0.1,0.1); # prior to push mass onto resolved genotypes
        }
    }
    abundance ~ uniform(0,1); # this could be changed to something more clever
    
    for( i in 1:L ){
        vector[4] theta;
        theta <- rep_vector(0,4);
        for( k in 1:K ){
            vector[4] t_k;
            # probability that strain k has genotype 0,0
            t_k[1] <- (1-genotype[k,linksites[i,1]]) * (1-genotype[k,linksites[i,2]]);
            # same for genotype 0,1; then 1,0; then 1,1
            t_k[2] <- (1-genotype[k,linksites[i,1]]) * genotype[k,linksites[i,2]];
            t_k[3] <- genotype[k,linksites[i,1]] * (1-genotype[k,linksites[i,2]]);
            t_k[4] <- genotype[k,linksites[i,1]] * genotype[k,linksites[i,2]];
            theta <- theta + abundance[k] * t_k; # weight the probabilities by strain abundance
        }
        increment_log_prob( multinomial_log(linkcounts[i],theta) );
    }
}

"
sm <- stan_model(model_name="3C_genotypes",model_code=genotype3C)


#
# simulate some data to test model fit
#
V<-1000  # number of variant sites
K<-2     # number of strains
geno<-matrix(nrow=K,ncol=V)
abund <- runif(K)
#abund <- c(0.5,0.5) # diploid with two haplotypes

for(v in 1:V){
    # choose a zero strain
    zero <- ceiling(runif(1,0,K))
    # choose a one strain from the rest
    one <- ceiling(runif(1,0,K-1))
    if(one >= zero){ one <- one + 1 }
    geno[zero,v] <- 0
    geno[one,v] <- 1
    # set the rest randomly
    for(k in 1:K){
        if(k == zero || k == one){ next }
        geno[k,v] <- rbinom(1,1,0.5)
    }
}

# simulate the links
readpairs=100000
ligation_error_rate=0.01
unif_weight=0.25
exp_mean=20
allpairs <- array(rep(0,V*V*2*2),dim=c(V,V,2,2))
for(r in 1:readpairs){
    site_i <- 0
    site_j <- 0
    geno_i <- -1
    geno_j <- -1
    if(runif(1)<ligation_error_rate){
        # choose two sites to link as a result of a proximity ligation error
        site_i <- ceiling(runif(1,0,V))
        site_j <- ceiling(runif(1,0,V))
        # choose two strains to link
        strains <- rmultinom(1,size=2,prob=abund)
        strain1 <- which(strains!=0)[0]
        strain2 <- tail(which(strains!=0),n=1)
        # add a link
        allpairs[site_i,site_j,geno[strain1,site_i]+1,geno[strain2,site_j]+1] = allpairs[site_i,site_j,geno[strain1,site_i]+1,geno[strain2,site_j]+1] + 1
    }else{
        # choose a strain to work with
        strain <- rmultinom(1,size=1,prob=abund)
        strain <- which(strain!=0)
        # choose the first site
        site_i <- ceiling(runif(1,0,V))
        # choose the 2nd site at a distance given by an exponential+uniform mixture
        dist_j = 0
        if(runif(1)<unif_weight){
            dist_j <- ceiling(runif(1,0,V))
        }else{
            dist_j <- ceiling(rexp(1,1/exp_mean))
        }
        site_j = (site_i + dist_j - 1) %% V + 1
        # link the sites
        allpairs[site_i,site_j,geno[strain,site_i]+1,geno[strain,site_j]+1] = allpairs[site_i,site_j,geno[strain,site_i]+1,geno[strain,site_j]+1] + 1
    }
}

linkcounts<-matrix(nrow=sum(allpairs>0),ncol=4)
linksites<-matrix(nrow=sum(allpairs>0),ncol=2)
c <- 1
for(i in 1:V){
    for(j in i:V){
        newrow <- c(allpairs[i,j,1,1],allpairs[i,j,1,2],allpairs[i,j,2,1],allpairs[i,j,2,2])
        newrow2 <- c(allpairs[j,i,1,1],allpairs[j,i,1,2],allpairs[j,i,2,1],allpairs[j,i,2,2])
        newrow <- newrow + newrow2
        if(sum(newrow)==0){ next }
        linkcounts[c,] <- newrow
        linksites[c,] <- c(i,j)
        c <- c+1
    }
}
linkcounts <- linkcounts[1:c-1,]
linksites <- linksites[1:c-1,]
hist(linksites[linksites[,1]==1,2],breaks=100) # to see if the link distance distribution looks reasonable

#
# END test data simulation
#


# fit the model to data
dat <- list(V=V, L=nrow(linkcounts), linkcounts=linkcounts, linksites=linksites, K=K);
vbfit <- vb(sm, data = dat, tol_rel_obj=0.001, iter=5000, sample_file = 'geno.csv', algorithm="meanfield")


# find the row permutation with minimum RMSE and return it
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

# measure model accuracy
infgeno<-t(matrix(get_posterior_mean(vbfit,pars="genotype"),V,K))
infgeno <- bestorder(infgeno, geno)
rinfgeno <- bestorder(round(infgeno), geno) # measure RMSE of the MAP estimate

blah <- matrix(rbinom(K*V,1,0.5),nrow=K,ncol=V)
blah <- bestorder(blah, geno)

