#9/18/19, by Doc Edge
#Goal: Take in a set of IBD calls using a "probe" dataset and the real one on which
#it's based, analyzing the IBD between probe and non-probe datasets. Uses GenomicRanges.


library(GenomicRanges)

#1.9 cM region
dat1 <- read.table("IBD_real_and_probe.chr19.lb45136023.ub45817621.vcf.gz.ibd", as.is = TRUE) 
#ids <- unique(c(dat1$V1, dat1$V3))

#Keep only matches where one member is a probe and the other is not
V1.is.probe <- grepl(".probe", dat1[,1])
V3.is.probe <- grepl(".probe", dat1[,3])
dat1 <- dat1[V1.is.probe + V3.is.probe == 1,]

#Remove matches where the probe is built from the same genome.
V1.orig.ID <- gsub(".probe", "", dat1[,1])
V3.orig.ID <- gsub(".probe", "", dat1[,3])
dat1 <- dat1[V1.orig.ID != V3.orig.ID,]

#Finally, remove according to length and LOD threshold
#dat1 <- dat1[dat1[,8] > 3,] #LOD score cutoff
dat1 <- dat1[dat1[,9] > 1,] #cM cutoff

#How many of these include the target site?
mean(dat1[,6] < 45411941 & dat1[,7] > 45411941)

mean(dat1[,6] < 45411941)
mean(dat1[,7] > 45411941)


#IBS.ranges<-makeGRangesFromDataFrame(dat,seqnames.field = "V5",start.field="V6",end.field="V7",ignore.strand=TRUE, keep.extra.columns = TRUE)



#5.9 cM region
dat3 <- read.table("IBD_real_and_probe.chr19.lb43834377.ub47047801.vcf.gz.ibd", as.is = TRUE) 
#ids <- unique(c(dat3$V1, dat3$V3))

#Keep only matches where one member is a probe and the other is not
V1.is.probe <- grepl(".probe", dat3[,1])
V3.is.probe <- grepl(".probe", dat3[,3])
dat3 <- dat3[V1.is.probe + V3.is.probe == 1,]

#Remove matches where the probe is built from the same genome.
V1.orig.ID <- gsub(".probe", "", dat3[,1])
V3.orig.ID <- gsub(".probe", "", dat3[,3])
dat3 <- dat3[V1.orig.ID != V3.orig.ID,]

#Finally, remove according to length and LOD threshold
#dat3 <- dat3[dat3[,8] > 3,] #LOD score cutoff
dat3 <- dat3[dat3[,9] > 3,] #cM cutoff

#How many of these include the target site?
mean(dat3[,6] < 45411941 & dat3[,7] > 45411941)

mean(dat3[,6] < 45411941)
mean(dat3[,7] > 45411941)


unique(dat1[,1])
unique(dat3[,1])

haps.1 <- paste(dat1[,1], dat1[,2], sep = ".")
haps.3 <- paste(dat3[,1], dat3[,2], sep = ".")
length(unique(haps.1[dat1[,6] < 45411941 & dat1[,7] > 45411941]))
length(unique(haps.3[dat3[,6] < 45411941 & dat3[,7] > 45411941]))

IBS.1<-makeGRangesFromDataFrame(dat1,seqnames.field = "V5",start.field="V6",end.field="V7",ignore.strand=TRUE, keep.extra.columns = TRUE)

IBS.3<-makeGRangesFromDataFrame(dat3,seqnames.field = "V5",start.field="V6",end.field="V7",ignore.strand=TRUE, keep.extra.columns = TRUE)


#plan to make plot: for each haplotype covered by at least one probe,
#produce a reduced set of ranges. Keep track of the starts and ends of these ranges.
#to make the plot, just marge along and see how many starts are to the left,
#then subtract off the number of ends tthat theo the left.
#this function will rely on the fact that the non-probe match is always in V1.
starts.and.ends <- function(IBSranges){
	haps <- paste(IBSranges$V1, IBSranges$V2, sep = ".")
	unique.haps <- unique(haps)
	starts <- numeric(0)
	ends <- numeric(0)
	for(i in 1:length(unique.haps)){
		hap <- unique.haps[i]
		ranges.hap <- reduce(IBSranges[haps == hap,])
		starts <- c(starts, start(ranges.hap))
		ends <- c(ends, end(ranges.hap))
	}
	cbind(starts, ends)
}

start.end.1 <- starts.and.ends(IBS.1)
start.end.3 <- starts.and.ends(IBS.3)

location <- seq(253938, 59097160, length.out = 10000)

covered.1 <- sapply(location, FUN = function(x, s.e){sum(x > s.e[,1]) - sum(x > s.e[,2]) }, s.e = start.end.1 )

covered.3 <- sapply(location, FUN = function(x, s.e){sum(x > s.e[,1]) - sum(x > s.e[,2]) }, s.e = start.end.3 )

pal <- c('#e66101','#fdb863','#b2abd2','#5e3c99')

pdf("probe.pdf", width = 5, height = 3)
par(las = 1, mar = c(3.2, 3.2, .6, .6), mgp = c(2.2, .7, 0), cex.lab = .9, cex.axis = .9)
plot(location, covered.1/(872*2), type = "l", bty = "n", ylim = c(0, 1), xlim = c(45411941 - 5e6, 45411941 + 5e6), lwd = 1.5, ylab = "Proportion of haplotypes covered", xlab = "Chromosomal location (mB)", xaxt = "n", col = pal[1], yaxt = "n")
axis(1, at = 45411941 + c(-5, -3, -1, 1, 3, 5)*1e6, labels = c("-5", "-3", "-1", "1", "3", "5"))
axis(2, at = c(0, .2, .4, .6, .8, 1), labels = c("0", ".2", ".4", ".6", ".8", "1"))
lines(location, covered.3/(872*2), col = pal[2], lwd = 1.5)
lines(rep(45136023,2), c(0,1), lty = 2, lwd = .5, col = pal[1])
lines(rep(45817621,2), c(0,1), lty = 2, lwd = .5, col = pal[1])
lines(rep(43834377,2), c(0,1), lty = 2, lwd = .5, col = pal[2])
lines(rep(47047801,2), c(0,1), lty = 2, lwd = .5, col = pal[2])
legend("topright", legend = c(1, 3), col = pal[1:2], lty = 1, lwd = 1.5, bty = "n", title = "cM cutoff")
dev.off()


