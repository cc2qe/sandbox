# get the per segment min fraction distribution

# pwd
# # /mnt/thor_pool1/user_data/cc2qe/projects/cin/e1_2013-09-27
# cat TCGA-*/*.subpops.hybrid.txt | grep -v "^#" | cut -f 12 > obs_min_frac.txt

mf <- scan('/Users/cc2qe/Documents/hall_lab/cin/e1_2013-09-27/obs_min_frac.txt', what='raw')
class(mf) <- 'numeric'

hist(mf, breaks=80, col='steelblue3', prob=TRUE, axes=FALSE, main='Per segment minor allele fraction of tumors')
d <- density(mf)
d$y.max
peak.x <- d$x[d$y == max(d$y)] # this returns the x which has the maximum value of y.
lines(d, col='red', lwd=2)
axis(2)
axis(1, at=seq(0,0.5,0.05))

abline(v=peak.x)

t <- rnorm(100)
thist <- hist(t, plot=TRUE)

barplot(thist$density, col='steelblue3', space=0, horiz=TRUE)
axis(2)
axis(2, at=thist$mids)