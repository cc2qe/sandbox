# Draw a histogram from a text file
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
x <- as.numeric(scan(file))

pdf('histogram.pdf', height=4, width=5)
hist(x, col='steelblue3', breaks=50, main=paste0('Histogram of ', args[1]), xlab=args[1])
dev.off()