#!/usr/bin/env Rscript
library(ape)
fff<-commandArgs(TRUE)[1]
ttt<-read.nexus(fff)
bbb<-balance(ttt)
# a balanced topology will have two clades with two descendants
write(paste(fff,sum(bbb==2)==2,sep="\t"),file="")

