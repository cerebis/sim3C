#!/usr/bin/env Rscript
library(ape)

ttt<-read.nexus(commandArgs(TRUE)[1])
bbb<-balance(ttt)
# a balanced topology will have two clades with two descendants
sum(bbb==2)==2

