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
        d = line[:-1].split("\t")
        num_strains = len(d)-1
        true_sites[d[0]] = d[1:]
    return true_sites

def parse_sites(sites_filename):
    global site_ids
    global num_sites
    global num_samples
    global num_strains
    sites_file = open(sites_filename)
    for line in sites_file:
        line = line[:-1]
        if line[0:3] == "U<-":
            num_sites = int(line[3:])
        if line[0:3] == "T<-":
            num_samples = int(line[3:])
        if line[0:3] == "S<-":
            num_strains = int(line[3:])
        if not line[0:7] == "siteids":
            continue
        line = line[13:]
        line = line[:-1]
        d = line.split(",")
        for i in range(len(d)):
            site_ids[i]=d[i]

    print "parsed " + str(len(site_ids)) + " site ids\n"

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
        return

    for i in range(num_strains):
        v = visited.copy()
        o = order[:]
        if i in v:
            continue
        v[i] = 1
        o.append(i)
        permuter(v, o)


def compute_accuracy( order ):
    t = -1
    dist = 0.0
    global best_dist
    global best_order
    global true_sites
    global site_ids
    unknown = 0
    total = 0
    for s in order:
        t += 1
        for j in range(num_sites):
            # get the truth for this site
            truth = [0,0,0,0]
            total += 1
            if not site_ids[j] in true_sites:
                unknown += 1
                continue
            truth[alpha_index[true_sites[site_ids[j]][t]]] = 1
            
            ssd = 0
            for i in range(len(truth)):
                ssd += pow(inferred[i][s][j]-truth[i], 2)
            dist += ssd

    dist = pow(dist, 0.5)
    if dist < best_dist:
        best_dist = dist
        best_order = order

#    print "unknown " + str(unknown)
#    print "total " + str(total)
#    print "dist " + str(dist)
#    exit(0) 

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Measure accuracy of a genotype reconstruction relative to a ground truth SNV matrix')
    parser.add_argument('--bpnmf', required=True, help='Path to bpnmf output')
    parser.add_argument('--truth', required=True, help='Path to true genotypes file')
    parser.add_argument('--sites', required=True, help='Path to R data file with site IDs')
    args = parser.parse_args()

    parse_sites(args.sites)
    inferred = parse_bpnmf(args.bpnmf)
    true_sites = parse_truth(args.truth)
    print "len ts " + str(len(true_sites))
    print "len inferred[0] " + str(len(inferred[0]))
    permuter(dict(), [])
    print "Best ordering found is " + ",".join(map(str,best_order))
    print "Best distance is " + str(best_dist)

    t = -1
    for s in best_order:
        t += 1
        for j in range(num_sites):
            if not site_ids[j] in true_sites:
                continue
            inf = ""
            for i in range(4):
                inf += ", " + str(inferred[i][s][j])
            print "strain " + str(t) + " site " + site_ids[j] + " true is " + true_sites[site_ids[j]][t] + " inferred " + inf


