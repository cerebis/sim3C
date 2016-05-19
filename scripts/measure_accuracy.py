#!/usr/bin/env python


import argparse
import os
import subprocess

best_dist = 99999999999999.0
best_order = []
num_strains = 0
num_samples = 0
true_sites = dict()
site_ids = dict()

alpha_index = {'A':0,'C':1,'G':2,'T':3}

def parse_truth(truth_filename):
    truth_file = open(truth_filename)
    for line in truth_file:
        d = line.split("\t")
        num_strains = len(d)-1
        true_sites[d[0]] = d[1:]
    return true_sites

def parse_bpnmf(bpnmf_filename):
    bpnmf_file = open(bpnmf_filename)
    global num_strains
    # init 2D array of tip partials
    tip_partials = [[[0 for x in range(num_sites)] for s in range(num_strains)] for x in range(len(alpha_index))]

    ll = -2 # skip the first line (csv header)
    for line in bpnmf_file:
        if line.startswith("#"):
            continue 
        ll += 1
        if ll < 0:
            continue
        d = line.split(",")
        for i in range(len(alpha_index)):
            for s in range(num_strains):
                begin = num_samples + 1 + num_sites * num_strains * i + num_sites * s
                end = begin + num_sites
                for j in range(begin,end):
                    tip_partials[i][s][j-begin] += float(d[j])

        # normalize to a tip partial distribution
        for s in range(num_strains):
            for j in range(num_sites):
                m = 0
                for i in range(len(alpha_index)):
                    m = m + tip_partials[i][s][j]
                for i in range(len(alpha_index)):
                    tip_partials[i][s][j] = tip_partials[i][s][j] / m

    return tip_partials

# calculate accuracy as minimum euclidean distance among all strain permutations

def permuter( visited, order ):
    if(len(visited)==num_strains):
        compute_accuracy(order)

    v = visited.copy()
    o = order[:]
    for i in range(num_strains):
        if i in v:
            continue
        v[i] = 1
        o.append(i)
        permuter(v, o)


def compute_accuracy( order ):
    t = 0
    dist = 0.0
    global best_dist
    global best_order
    global true_sites
    global site_ids
    for s in order:
        for j in range(num_sites):
            # get the truth for this site
            truth = [0,0,0,0]
            truth[alpha_index[true_sites[site_ids[j]][s]]] = 1
            
            ssd = 0
            for i in range(len(truth)):
                ssd += pow(inferred[s][j][i]-truth[s][j][i], 2)
            dist += ssd

    dist = pow(dist, 0.5)
    if dist < best_dist:
        best_dist = dist
        best_order = order


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Measure accuracy of a genotype reconstruction relative to a ground truth SNV matrix')
    parser.add_argument('--bpnmf', required=True, help='Path to bpnmf output')
    parser.add_argument('--truth', required=True, help='Path to true genotypes file')
    parser.add_argument('--num-sites', required=True, help='Number of variant sites found')
    parser.add_argument('--num-strains', required=True, help='Number of strains inferred')
    parser.add_argument('--num-samples', required=True, help='Number of samples taken')
    args = parser.parse_args()

    num_sites = int(args.num_sites)
    num_strains = int(args.num_strains)
    num_samples = int(args.num_samples)
    inferred = parse_bpnmf(args.bpnmf)
    true_sites = parse_truth(args.truth)
    print "len ts " + str(len(true_sites))
    print "len inferred[0] " + str(len(inferred[0]))
    permuter(dict(), [])
    print "Best ordering found is " + ",".join(best_order) + "\n"
    print "Best distance is " + str(best_dist)


