setwd('/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27')

# args <- commandArgs(trailingOnly=TRUE)
file <- '/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27/TCGA-A1-A0SJ/TCGA-A1-A0SJ-01.subpops.hybrid.ext.txt'
# file <- args[1]
outroot <- '/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27/test'
# outroot <- args[2]

s <- matrix(scan(file, what='raw', sep='\t'), ncol=13, byrow=TRUE)
rownames(s) <- s[,1]
colnames(s) <- s[1,]
s <- s[2:nrow(s), 2:ncol(s)]
class(s) <- 'numeric'

# defines the sizes of each chromosome, so that we can create running sums
chromSizes <- c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
names(chromSizes) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

chromRunsum <- c(0)
for (i in 1:length(chromSizes)) {
	chromRunsum[i+1] <- sum(chromSizes[1:i])
}
names(chromRunsum) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "ALL")


# plot the histogram of mutant population fractions
# make a vector of the subpopulation frequencies
q <- vector()
for (i in 1:nrow(s)) {
	# length <- s[i,2] - s[i,1]
	# print(length)
	
	q <- c(q, s[i,8])
}

# v is the version that accounts for the length of each segment, protecting from instances where the segmenter shreds a region of the genome.
v <- vector()
for (i in 1:nrow(s)) {
	length <- s[i,2] - s[i,1]
	k <- length %/% 100000
	# print(k)
	v <- c(v, rep(s[i,8], k))
}

# -------------------------------------------
# sandbox

library('mixtools')

match.model <- function(q, numClust) {
	qhist <- hist(q, breaks=seq(0,1,0.001), col='steelblue3', prob=TRUE, las=1, plot=FALSE)

	if (numClust < 2) {
		mixmdl <- NULL
		mixmdl$mu <- mean(v)
		mixmdl$sigma <- sd(v)
		mixmdl$lambda <- 1
		

	} else {
		mixmdl <- normalmixEM(q, k=numClust)
		# plot(mixmdl, which=2)
	}
	
	y <- NULL
	for (i in 1:numClust) {
		res <- 10000
		l <- mixmdl$lambda[i]
		mu <- mixmdl$mu[i]
		s <- mixmdl$sigma[i]
		
		y <- c(y, rnorm(res * l, mean=mu, sd=s))
	}
	
	# y <- y[y>=0 & y<=1]
	# min(y)
	yhist <- hist(y, breaks=seq(min(y)-0.001,max(y)+0.001,0.001), col='steelblue3', prob=TRUE, las=1, plot=FALSE)
	

	
	
	# yhist$mids
	# dev.new()
	# plot(qhist$mids, z, type='l', col='navy', main='cdf and ks distance', axes=FALSE)
	# lines(qhist$mids, p, type='l', col='red')
	# axis(1)
	# axis(2, at=seq(0,1,0.1), las=1)
	
	output <- NULL
	output$ks <- ks.test(v, y)
	output$mixmdl <- mixmdl
	output$qhist <- qhist
	output$yhist <- yhist
	
	return(output)
}

# make a function out of the quantile, so we can access the same data points in the test and model cdfs.
get.cd <- function(x, y, quantile) {
	if (length(x) != length(y)) {
		print("error: length of x and y not equal")
		return()
	}
	ind <- max(which(x<=quantile, arr.ind=TRUE))
	return(y[ind])
}

# riemann sum of the difference between curves
riemann.sum <- function(a.x, a.y, b.x, b.y, step) {
	# function is only defined where there is both a and b curves
	global.start <- max(min(a.x), min(b.x)) 
	global.end <- min(max(a.x), max(b.x))
	
	# the riemann sum
	r.sum <- 0
	for (mid in seq(global.start, global.end, step)) {
		diff <- abs(get.cd(a.x, a.y, mid) - get.cd(b.x, b.y, mid))
		r.sum <- r.sum + step * diff
	}
	return(r.sum)
}

count.clusters <- function(M) {
	bestClust <- 0
	
	# the differential between each successive cluster
	d <- NULL
	for (i in 2:length(M)) {
		d[i] <- M[i-1] - M[i]
	}
	
	# the total difference between the first and last point
	total <- sum(d[2:length(d)])
	
	# get the last element that accounts for at least 30% of the total drop OR the element of max drop
	for (i in length(d):2) {
		if (d[i] > 0.3 * total) {
			bestClust <- i
			break
		}
	}
	if (bestClust == 0) {
		bestClust <- rev(order(d))[1]
	}
	
	return(bestClust)
}

# calculate the cin score
get.cinScore <- function(x, mixmdl) {
	cin.score <- 0
	numVars <- length(x)
	for (i in 1:length(mixmdl$lambda)) {
		cin.score <- cin.score + mixmdl$lambda[i] * mixmdl$sigma[i]
	}
	cin.score <- cin.score * numVars
	
	return(cin.score)
}

# flatness score
flatness.sum <- function(a.x, a.y, step) {
	# function is only defined where there is both a and b curves
	global.start <- min(a.x)
	global.end <- max(a.x)
	
	# the riemann sum. (the other curve here is just a diagonal line y=x
	f.sum <- 0
	for (mid in seq(global.start, global.end, step)) {
		diff <- abs(get.cd(a.x, a.y, mid) - get.cd(a.x, a.x, mid))
		f.sum <- f.sum + step * diff
	}
	return(f.sum)
}

# --------------------------------------------
# calculate flatness

qhist <- hist(q, breaks=seq(0,1,0.001), col='steelblue3', prob=TRUE, las=1, plot=FALSE)
z <- NULL
for (i in 1:length(qhist$density)) {
	z[i] <- sum(qhist$density[1:i])/1000
}
flatness <- 1 - (2 * flatness.sum(qhist$mids, z, 0.001))

f.res <- 100
rhist <- hist(q, breaks=seq(0,1,1/f.res), col='steelblue3', prob=TRUE, las=1, plot=FALSE)
z.sort <- NULL
dens.sort <- sort(rhist$density, decreasing=TRUE)
for (i in 1:length(dens.sort)) {
	z.sort[i] <- sum(dens.sort[1:i])/f.res
}
flatness.sort <- 1 - (2 * flatness.sum(rhist$mids, z.sort, 1/f.res))

# plot flatness
pdf(paste0(outroot,'_flatness.pdf'), width=8, height=8)
plot(qhist$mids, rep(1,length(qhist$mids)), type='h', col='steelblue3', lwd=1.2, xlim=c(0,1), ylim=c(0,1), las=1, main=paste0('flatness = ', round(flatness,4)), ylab='cumulative density', xlab='quantile')
lines(qhist$mids, z, type='h', lwd=2, col='white')
lines(qhist$mids, z, type='l', lwd=2, col='black')
abline(a=0, b=1, col='red', lwd=2)
dev.off()

# plot sorted flatness
pdf(paste0(outroot,'_flatness.sort.pdf'), width=8, height=8)
plot(rhist$mids, rep(1,length(rhist$mids)), type='h', col='steelblue3', lwd=4, xlim=c(0,1), ylim=c(0,1), las=1, main=paste0('flatness.sort = ', round(flatness.sort,4)), ylab='cumulative density', xlab='quantile')
lines(rhist$mids, z.sort, type='h', lwd=5, col='white')
lines(rhist$mids, z.sort, type='l', lwd=2, col='black')
abline(a=0, b=1, col='red', lwd=2)
dev.off()

# write out flatness score
write(flatness, paste0(outroot,'_flatnessscore.txt'), ncol=1)

# write out flatness.sort score
write(flatness, paste0(outroot,'_flatnessSortScore.txt'), ncol=1)








