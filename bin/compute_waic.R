#!/usr/bin/env Rscript
#
# meta-sweeper - for performing parametric sweeps of simulated
# metagenomic sequencing experiments.
# Copyright (C) 2016 "Matthew Z DeMaere"
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

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

write("#waic, standard error", file="")
write(paste(lpd_hat - phat_waic, waic_se, sep=","), file="")
