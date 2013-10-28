setwd('/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27')

# args <- commandArgs(trailingOnly=TRUE)
file <- '/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27/TCGA-A1-A0SJ/TCGA-A1-A0SJ-01.subpops.hybrid.ext.txt'
outfile <- '/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27/test.pdf'

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

# dev.new(width=16, height=11)
pdf(outfile, 18,10)
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

v <- vector()
for (i in 1:nrow(s)) {
	length <- s[i,2] - s[i,1]
	k <- length %/% 100000
	# print(k)
	v <- c(v, rep(s[i,8], k))
}

par(fig=c(0.8,1,0,0.25), new=TRUE)
# qhist <- hist(q, breaks=seq(0,1,0.01), plot=FALSE)
# barplot(qhist$density, col='navy', border='navy', space=0, horiz=TRUE, xlab='frequency')
# # axis(2, at=axTicks(2), labels=qhist$breaks[axTicks(2)+1], las=1)

vhist <- hist(v, breaks=seq(0,1,0.01), plot=FALSE)
barplot(vhist$density, col='navy', border='navy', space=0, horiz=TRUE, xlab='frequency')
# axis(2, at=axTicks(2), labels=qhist$breaks[axTicks(2)+1], las=1)

dev.off()

# -------------------------------------------
# sandbox

library('mixtools')

for (k in 2:5) {
	
	mixmdl <- normalmixEM(v, k=k)
	plot(mixmdl, which=2)
	
	seq(0,1,0.01)
	vhist <- hist(v, breaks=seq(0,1,0.01), col='steelblue3', prob=TRUE, las=1)
	
	# for (i in 1:length(mixmdl$mu)) {
		# x <- seq(min(v), max(v), length=100)
		# y <- mixmdl$lambda[i] * dnorm(x, mean=(mixmdl$mu[i]), sd=mixmdl$sigma[i])
		# lines(x, y, col=rainbow(length(mixmdl$mu))[i], lwd=2)
	# }
	
	y <- rep(0, 100)
	for (i in 1:length(mixmdl$mu)) {
		x <- seq(min(v), max(v), length=100)
		y <- y + (mixmdl$lambda[i] * dnorm(x, mean=(mixmdl$mu[i]), sd=mixmdl$sigma[i]))
	}
	# y <- y/length(mixmdl$mu)
	
	lines(x, y, type='l', lwd=2, col='red')
	ks <- ks.test(vhist$density, y)
	print(ks)
}

plot(1, type='n', ylim=c(.8,1), xlim=c(10,16))
lines(ecdf(vhist$density), col='red', type='l')
lines(ecdf(y), col='blue', type='l')

# plot(vhist$density, type='h', col='navy', lwd=2)

maxDev <- 0
for (i in 1:length(vhist$density)) {
	if (abs(vhist$density[i] - y[i]) > maxDev) {
		maxDev <- abs(vhist$density[i] - y[i])
	}
}

maxDev
a <- abs(vhist$density - y)==maxDev

vhist$density[15]
y[15]

mixmdl$sigma
names(k)

# -------------------------------------------

match.model <- function(v, numClust) {
	vhist <- hist(v, breaks=seq(0,1,0.001), col='steelblue3', prob=TRUE, las=1, plot=FALSE)

	if (numClust < 2) {
		mixmdl <- NULL
		mixmdl$mu <- mean(v)
		mixmdl$sigma <- sd(v)
		mixmdl$lambda <- 1
		
		hist(v, prob=TRUE, main='Density Curves')
		xfit <- seq(min(v), max(v), length=1000)
		yfit <- dnorm(xfit, mean=mixmdl$mu, sd=mixmdl$sigma)
		lines(xfit, yfit, col='red', lwd=2)
	} else {
		mixmdl <- normalmixEM(v, k=numClust)
		plot(mixmdl, which=2)
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
	
	z <- NULL
	for (i in 1:length(vhist$density)) {
		z[i] <- sum(vhist$density[1:i])/1000
	}
	p <- NULL
	for (i in 1:length(yhist$density)) {
		p[i] <- sum(yhist$density[1:i])/1000
	}
	
	
	yhist$mids
	dev.new()
	plot(vhist$mids, z, type='l', col='navy', main='cdf and ks distance', axes=FALSE)
	lines(yhist$mids, p, type='l', col='red')
	axis(1)
	axis(2, at=seq(0,1,0.1), las=1)
	
	output <- NULL
	output$ks <- ks.test(v, y)
	output$mixmdl <- mixmdl
	
	return(output)
}

# --------------------------------------------

D <- NULL
for (i in 1:8) {
	m <- match.model(v, i)
	D[i] <- m$ks$statistic
}
plot(D, type='l', col='navy', lwd=2)

# ------------------------------------------

plot(vhist$mids, yhist$mids)
for (s in seq(min(c(vhist$mids, yhist$mids)),max(c(vhist$mids, yhist$mids)),0.001)) {
	for (i in 1:length(vhist$mids)) {
		if (s < vhist$mids[
	}
}

maxDev <- 0
for (i in 1:length(z)) {
	if (abs(z[i]-p[i]) > maxDev) {
		maxDev <- abs(z[i]-p[i])
	}
}
maxDev
m <- abs(z-p)
lines(m, type='l', col='magenta', main='ks distance', ylab='abs(model - test)')
i <- which(abs(z-p)==maxDev, arr.ind=TRUE)
points(i, abs(z[i]-p[i]), pch=4)

length(yhist$density)

ks.test(v,y)

# plot(mixmdl, which=2)

# plot(p, typ='l')






lines(seq(0,1,length=100), vhist$density, type='l', col='red')
length(vhist$density)

y <- c(rnorm(10000*mixmdl$lambda[1], mean=mixmdl$mu[1], sd=mixmdl$sigma[1]), rnorm(1000*mixmdl$lambda[2], mean=mixmdl$mu[2], sd=mixmdl$sigma[2]))
ks <- ks.test(v, y)
ks

plot(y)

hist(y, breaks=600)

# --------------------------------------------






a <- ecdf(vhist$density)
b <- ecdf(y)
x <- seq(0,100,length=10000)
plot(x, a(x), type='l', col='red', xlim=c(0,20))
lines(x, b(x), type='l', col='blue')




ks <- ks.test(vhist$density, y) # this is wrong because i'm giving it pdf. I should give it the raw data (like v)
# ks <- ks.test(v, pnorm(seq(0,1,0.01), mean=0.1192611, sd=0.05766979))
ks <- ks.test(v, c(rnorm(1000*mixmdl$lambda[1], mean=mixmdl$mu[1], sd=mixmdl$sigma[1]), rnorm(1000*mixmdl$lambda[2], mean=mixmdl$mu[2], sd=mixmdl$sigma[2])))
ks

?ks.test
a(3)
b(3)

# ------------------------------------
# ks example

require(graphics)

x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)
# Does x come from a shifted gamma distribution with shape 3 and rate 2?
ks.test(x+2, "pgamma", 3, 2) # two-sided, exact
ks.test(x+2, "pgamma", 3, 2, exact = FALSE)
ks.test(x+2, "pgamma", 3, 2, alternative = "gr")

# test if x is stochastically larger than x2
x2 <- rnorm(50, -1)
plot(ecdf(x), xlim = range(c(x, x2)))
plot(ecdf(x2), add = TRUE, lty = "dashed")
t.test(x, x2, alternative = "g")
wilcox.test(x, x2, alternative = "g")
ks.test(x, x2, alternative = "l")



