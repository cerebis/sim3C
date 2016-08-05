#!/usr/bin/env Rscript
stanfile <- commandArgs(TRUE)[1]
U <- as.numeric(commandArgs(TRUE)[2])
T <- as.numeric(commandArgs(TRUE)[3])

abab<-read.table(stanfile,sep=",",header=TRUE)
ff <- ncol(abab) - U*T + 1
ffl <- ncol(abab)

lpd_hat<-sum(log(colMeans(exp(abab[,ff:ffl]))))
phat_waic <- sum(apply(abab[,ff:ffl],2,var))


elpd_hats <- log(colMeans(exp(abab[,ff:ffl]))) - apply(abab[,ff:ffl],2,var)
waic_se <- sqrt(nrow(abab) * var(elpd_hats))

write("waic, standard error", file="")
write(paste(lpd_hat - phat_waic, waic_se, sep=","), file="")
