#!/usr/bin/env Rscript

# print usage
usage <- function() {
	cat(
'usage: z_prob.R <z> <mean> <stdev>
z_prob.R
author: Colby Chiang (cchiang@wusm.wustl.edu)
description: calculate the probability of a z-score
positional arguments:
  z        z-score to test
  mean     mean of normal distribution
  stdev    standard deviation of normal distribution
')
}

args <- commandArgs(trailingOnly=TRUE)

# Check input args
if (length(args) < 3) {
	usage()
	quit(save='no', status=1)
}

z <- as.numeric(args[1])
mean <- as.numeric(args[2])
stdev <- as.numeric(args[3])



cat(pnorm(z, mean=mean, sd=stdev), sep='\n')