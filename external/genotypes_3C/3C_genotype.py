#!/usr/bin/env python
import itertools
import pystan
import numpy as np
import argparse
from datetime import datetime



def rmse(pred, truth):
    return np.sqrt(((pred - truth) ** 2).mean())


def get_posteriors(fit, K, V):
    with open(fit['args']['sample_file'], 'r') as inh:
        header = None
        table = []
        for line in inh:
            if line.startswith('#') or line.startswith('Stepsize') or line.startswith('eta'):
                continue
            if not header:
                header = line.split(',')
                continue
            table.append([float(tok) for tok in line.split(',')])
        table = np.array(table)
        genos = np.ndarray(shape=(K, V))
        for k in xrange(K):
            genos[k] = table[:, k + 1:-K:K].mean(axis=0)
        return {'geno': genos,
                'abund': table[:, -K:].mean(axis=0)}


genotype3C = """
data {
    # number of sites (nodes in graph)
    int V;
    # number of links among sites (edges in graph)
    int L;
    # counts of SNV co-observation. Genotypes are binary variables: mutant or not. Four possible combinations among two sites: 00,01,10,11 
    int linkcounts[L, 4];
    # the indices of the two sites involved in each link
    int linksites[L, 2];
    # number of components to factorize
    int K; 
}
parameters {
    # the genotype of the K strains. Ideally this would be a binary (or base 4) variable
    real<lower=0, upper=1> genotype[K, V];
    
    # the abundance of each strain
    simplex[K] abundance;
}
model {
    for(k in 1:K){
        for(v in 1:V){
            # prior to push mass onto resolved genotypes
            genotype[k,v] ~ beta(0.3, 0.3);
        }
    }
    
    # this could be changed to something more clever
    abundance ~ uniform(0, 1);
    
    for( i in 1:L ){

        vector[4] theta;
        theta = rep_vector(0, 4);

        for( k in 1:K ){

            vector[4] t_k;
            
            # probability that strain k has genotype 0,0
            t_k[1] = (1 - genotype[k, linksites[i, 1]]) * (1 - genotype[k, linksites[i, 2]]);
            
            # same for genotype 0,1; then 1,0; then 1,1
            t_k[2] = (1 - genotype[k, linksites[i, 1]]) *      genotype[k, linksites[i, 2]];
            t_k[3] =      genotype[k, linksites[i, 1]]  * (1 - genotype[k, linksites[i, 2]]);
            t_k[4] =      genotype[k, linksites[i, 1]]  *      genotype[k, linksites[i, 2]];
            
            # weight the probabilities by strain abundance
            theta = theta + abundance[k] * t_k;
        }

        target += multinomial_lpmf(linkcounts[i] | theta);
    }
}
"""

parser = argparse.ArgumentParser(description='VB genotypes')
parser.add_argument('--num-sites', metavar='INT', required=True, type=int, help='Number of SVN sites')
parser.add_argument('--num-geno', metavar='INT', required=True, type=int, help='Number of genotypes')
parser.add_argument('--num-pairs', metavar='INT', required=True, type=int, help='Number of simulated ligation read-pairs')
parser.add_argument('--error-rate', metavar='FLOAT', type=float, default=0.01, help='Rate of spurious ligation products [0.01]')
parser.add_argument('--seed', metavar='INT', type=int, default=datetime.now().microsecond, help='Random seed')
parser.add_argument('output', help='output file')
args = parser.parse_args()

V = args.num_sites
K = args.num_geno
read_pairs = args.num_pairs
ligation_error_rate = args.error_rate

UNIF_WEIGHT = 0.25
EXP_MEAN = 20.0

rnd_state = np.random.RandomState(args.seed)


geno = np.zeros((K, V), dtype=int)
abund = rnd_state.uniform(size=K)
abund /= abund.sum()

# create genotypes
elems = set(range(K))
for v in xrange(V):
    xi = rnd_state.choice(K, replace=False, size=2)
    geno[xi, v] = [0, 1]
    if K > 2:
        xi = rnd_state.choice(list(elems - set(xi)), size=K - 2)
        geno[xi, v] = rnd_state.binomial(1, 0.5, size=K - 2)

# create linked sites
all_pairs = np.zeros((V, V, 2, 2), dtype=int)

for ri in xrange(read_pairs):
    if rnd_state.uniform() < ligation_error_rate:
        # spurious product
        site_i, site_j = rnd_state.randint(0, V, size=2)
        strain1, strain2 = rnd_state.choice(K, size=2, p=abund, replace=False)
        all_pairs[site_i, site_j, geno[strain1, site_i], geno[strain2, site_j]] += 1
    else:
        # cis product
        strain = rnd_state.choice(K, p=abund)
        site_i = rnd_state.randint(0, V)
        dist_j = None
        if rnd_state.uniform() < UNIF_WEIGHT:
            dist_j = rnd_state.randint(0, V)
        else:
            dist_j = np.ceil(rnd_state.exponential(EXP_MEAN)).astype(int)
        site_j = (site_i + dist_j - 1) % V
        all_pairs[site_i, site_j, geno[strain, site_i], geno[strain, site_j]] += 1

num_pairs = np.sum(all_pairs > 0)
print 'Number of connecting pairs', num_pairs

link_counts = np.zeros((num_pairs, 4), dtype=int)
link_sites = np.zeros((num_pairs, 2), dtype=int)

n = 0
for i in xrange(V):
    for j in xrange(i, V):
        new_row = all_pairs[i, j, :, :] + all_pairs[j, i, :, :]
        if new_row.sum() == 0:
            # not connecting
            continue
        link_counts[n, :] = new_row.flatten()
        link_sites[n, :] = (i, j)
        n += 1

# trim to initialised length
link_counts = link_counts[:n, ]
link_sites = link_sites[:n, ]
link_sites += 1

# plot histogram of link sites
# import matplotlib.pyplot as plt
# plt.hist(link_sites[link_sites[:,0]==1,1], bins=100)
# plt.show()

# run stan
dat = {'V': V, 'L': link_counts.shape[0], 'linkcounts': link_counts, 'linksites': link_sites, 'K': K}
sm = pystan.StanModel(model_name='a3C_genotypes', model_code=genotype3C)
vb_fit = sm.vb(data=dat, tol_rel_obj=0.001, iter=5000, sample_file='geno.csv',
               seed=args.seed, algorithm="meanfield", verbose=True)

infer = get_posteriors(vb_fit, K, V)
np.savetxt(args.output, np.hstack((infer['geno'].T, geno.T)), fmt='%.4e')

scores = {}
for ix in itertools.permutations(range(K)):
    scores[ix] = rmse(infer['geno'][np.array(ix)], geno)

print 'permuted scores against truth'
scores = sorted(scores.items(), key=lambda x: x[1])
for sc in scores:
    print '{0} rmse= {1:.5f}'.format(sc[0], sc[1])

print 'reord abund:', infer['abund'][np.array(scores[0][0])]
print 'truth abund:', abund

