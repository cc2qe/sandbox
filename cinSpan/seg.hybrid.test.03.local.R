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


# par(mfrow=c(4,1), mar=c(2.5,2,1,0.5), oma=c(0.5,2,.5,0.5))

# dev.new(width=18, height=10)
pdf(paste0(outroot,'.pdf'), 18,10)
par(mar=c(2,1.7,1,0), oma=c(0.5,2,.5,0.5))

par(fig=c(0,0.8,0.75,1), new=TRUE)
# plot segments (copy number)
plot(1, type='n', ylim=c(-2,2), main='log2 copy number (obs: blue, calc: green)', xlab='genome position', ylab='ylab', xlim=c(0,chromRunsum['ALL']), las=1)
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], log(s[i,9],2), chromRunsum[rownames(s)[i]]+s[i,2], log(s[i,9],2), lwd=3, col='navy')
}
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], log(s[i,10],2), chromRunsum[rownames(s)[i]]+s[i,2], log(s[i,10],2), lwd=3, col='seagreen3')
}
abline(h=0)
for (b in chromRunsum) {
	abline(v=b, col='gray50')
}

par(fig=c(0,0.8,0.5,0.75), new=TRUE)
# plot minor allele fraction
plot(1, type='n', ylim=c(0,0.5), main='minor allele fraction (obs: blue, calc: green)', xlab='genome position', ylab='ylab', xlim=c(0,chromRunsum['ALL']), las=1)
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], s[i,11], chromRunsum[rownames(s)[i]]+s[i,2], s[i,11], lwd=3, col='navy')
}
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], s[i,12], chromRunsum[rownames(s)[i]]+s[i,2], s[i,12], lwd=3, col='seagreen3')
}
# abline(h=0)
for (b in chromRunsum) {
	abline(v=b, col='gray50')
}

par(fig=c(0,0.8,0.25,0.5), new=TRUE)
# plot number of alleles in mutant
plot(1, type='n', ylim=c(0,6), main='variant num alleles (calc) (blue: total, red: minor)', xlab='genome position', ylab='ylab', xlim=c(0,chromRunsum['ALL']), las=1)
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], s[i,4], chromRunsum[rownames(s)[i]]+s[i,2], s[i,4], lwd=3, col='navy')
}
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], s[i,4]*s[i,6], chromRunsum[rownames(s)[i]]+s[i,2], s[i,4]*s[i,6], lwd=3, col='firebrick')
}
abline(h=2)
for (b in chromRunsum) {
	abline(v=b, col='gray50')
}

par(fig=c(0,0.8,0,.25), new=TRUE)
# plot subpop fraction
plot(1, type='n', ylim=c(0,1), main='variant pop fraction (calc)', xlab='genome position', ylab='ylab', xlim=c(0,chromRunsum['ALL']), las=1)
for (i in 1:nrow(s)) {
	segments(chromRunsum[rownames(s)[i]]+s[i,1], s[i,8], chromRunsum[rownames(s)[i]]+s[i,2], s[i,8], lwd=3, col='navy')
}
# abline(h=0)
for (b in chromRunsum) {
	abline(v=b, col='gray50')
}

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

par(fig=c(0.8,1,0,0.25), new=TRUE)
qhist <- hist(q, breaks=seq(0,1,0.01), plot=FALSE)
barplot(qhist$density, col='navy', border='navy', space=0, horiz=TRUE, xlab='frequency')
# axis(2, at=axTicks(2), labels=qhist$breaks[axTicks(2)+1], las=1)

# vhist <- hist(v, breaks=seq(0,1,0.01), plot=FALSE)
# barplot(vhist$density, col='navy', border='navy', space=0, horiz=TRUE, xlab='frequency')
# # axis(2, at=axTicks(2), labels=qhist$breaks[axTicks(2)+1], las=1)

dev.off()

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


# --------------------------------------------

D <- NULL
R <- NULL
for (j in 1:6) {
	m <- match.model(q, j)
	# plot(m$mixmdl, which=2)
	D[j] <- m$ks$statistic
	qhist <- m$qhist
	yhist <- m$yhist
	mixmdl <- m$mixmdl

	z <- NULL
	for (i in 1:length(qhist$density)) {
		z[i] <- sum(qhist$density[1:i])/1000
	}
	p <- NULL
	for (i in 1:length(yhist$density)) {
		p[i] <- sum(yhist$density[1:i])/1000
	}

	R[j] <- riemann.sum(qhist$mids, z, yhist$mids, p, 0.001)

	# plot the cdf overlay curves
	pdf(paste0(outroot,'cdf_k',j,'.pdf'), width=8, height=8)
	plot(qhist$mids, z, type='l', col='navy', main=paste0('cdf and ks distance (', j, ' clusters)'), axes=FALSE, xlab='quantile', ylab='cumulative density')
	lines(yhist$mids, p, type='l', col='red')
	axis(1)
	axis(2, at=seq(0,1,0.1), las=1)
	dev.off()
	
	# plot the histogram with overlaid density curve.
	pdf(paste0(outroot,'_hist_k',j,'.pdf'), width=8, height=8)
	xfit <- seq(min(q), max(q), length=1000)
	yfit <- rep(0, length(xfit))
	for (i in 1:j) {
		yfit <- yfit + mixmdl$lambda[i] * dnorm(xfit, mean=mixmdl$mu[i], sd=mixmdl$sigma[i])
	}
	hist(q, prob=TRUE, main=paste0('Density Curves (',j,' clusters)'), breaks=50, col='grey', las=1)
	lines(xfit, yfit, col='red', type='l', lwd=2)
	dev.off()

}

# output the best cluster
bestClust.ks <- count.clusters(D)
bestClust.riemann <- count.clusters(R)
# pdf(paste0(outroot,'_bestClust.pdf'), width=8, height=8)
# plot(1, type='n', col='navy', lwd=2, main=paste0('Best # of clusters: ', bestClust), xlab='', ylab='')
# dev.off()

# plot KS D statistic
pdf(paste0(outroot,'_ks_D.pdf'), width=8, height=8)
plot(D, type='l', col='navy', lwd=2, main='ks D statistic (red: ks opt, black: riemann opt)', xlab='# clusters', ylab='D')
points(bestClust.ks, D[bestClust.ks], pch=2, col='red')
points(bestClust.riemann, D[bestClust.riemann], pch=5, col='black')
dev.off()

# plot Riemann integral
pdf(paste0(outroot,'_riemann.pdf'), width=8, height=8)
plot(R, type='l', col='navy', lwd=2, main='Area between curves (red: ks opt, black: riemann opt)', xlab='# clusters', ylab='area between curves')
points(bestClust.ks, R[bestClust.ks], pch=5, col='red')
points(bestClust.riemann, R[bestClust.riemann], pch=5, col='black')
dev.off()


# now calculate it again with bestClust (based on riemann) and calculate the cin score
m <- match.model(q, bestClust.riemann)
D[j] <- m$ks$statistic
qhist <- m$qhist
yhist <- m$yhist
mixmdl <- m$mixmdl

cin.score <- get.cinScore(q, mixmdl)

write(cin.score, paste0(outroot,'_cinscore.txt'), ncol=1)









# -----------------------------------------
# xfit <- seq(min(q), max(q), length=1000)
# yfit <- rep(0, length(xfit))
# for (i in 1:3) {
	# yfit <- yfit + m$mixmdl$lambda[i] * dnorm(xfit, mean=m$mixmdl$mu[i], sd=m$mixmdl$sigma[i])
# }
# hist(v, prob=TRUE, main='Density Curves', breaks=50, col='grey', las=1)
# lines(xfit, yfit, col='red', type='l', lwd=2)
