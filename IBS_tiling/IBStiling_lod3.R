
#9/4/19 by Doc Edge
#Goal: Phase and call IBS to reproduce the IBS tiling results of Edge and Coop, "Attacks on genetic 
#privacy via uploads to genealogical databases." 
#This script reproduces the Appendix result that only includes segments with a LOD score for IBD > 3.
#This script requires:

#### R with the data.table and GenomicRanges packages installed





#########################################################



library(data.table)
library(GenomicRanges)


dat <- fread("ibd/IBD_HO_Euro_1kg_chr1.ibd") 
ids <- unique(c(dat$V1, dat$V3))
rm(dat)

#Strategy: make master GRanges for chrom,
#Subset into list w/ oneGRanges per individual
#Subset that list into haplotypes

#Loop through all the chromosomes, cM cutoffs, and bait numbers
chrs <- 1:22
cM.cuts <- c(1,3,5,8) 
bait.ns <- c(100, 200, 300, 400, 500, 600, 700, 800, 872)

set.seed(8675309)
bait <- sample(ids)

#Initialize arrays that store the tiling coverage for each person, chromosome at 
#every comparison sample size and cM reporting threshold
h1.covered <- array(dim = c(length(ids), length(chrs), length(bait.ns), length(cM.cuts)))
h2.covered <- h1.covered
IBD2.covered <- h1.covered

all.ranges.chr.list <- list()
rangelists.chr.list <- list()

#Reduces the GRanges where the other individual is in the bait set and the segments
#are longer than a cutoff cM.cut
reduce.in.bait <- function(GRangelist, bait, cM.cut){
	lapply(GRangelist, function(x){reduce(x[x$other.ind %in% bait & x$V9 > cM.cut])})
}


for(j in chrs){
	IBS<-read.table(paste("ibd/IBD_HO_Euro_1kg_chr", as.character(j), ".ibd", sep = ""),as.is=TRUE)

	IBS.ranges<-makeGRangesFromDataFrame(IBS,seqnames.field = "V5",start.field="V6",end.field="V7",ignore.strand=TRUE, keep.extra.columns = TRUE)

	all.ranges.chr.list[[j]] <- IBS.ranges

	IBS.rangelist <- GRangesList()
	for(i in 1:length(ids)){
		IBS.ind <- IBS.ranges[(IBS.ranges$V1 == ids[i] | IBS.ranges$V3 ==ids[i]) & IBS.ranges$V8 > 3,] #This is the only line that changes; filters out segments w/ LOD score less than 3.
		target.ind <- rep(ids[i], length(IBS.ind))
		IBS.ind$target.hap <- IBS.ind$V2
		IBS.ind$target.hap[IBS.ind$V3 == target.ind] <- IBS.ind$V4[IBS.ind$V3 == target.ind]
		IBS.ind$other.ind <- IBS.ind$V3
		IBS.ind$other.ind[IBS.ind$V3 == target.ind] <- IBS.ind$V1[IBS.ind$V3 == target.ind]
		IBS.ind$other.hap <- IBS.ind$V4
		IBS.ind$other.hap[IBS.ind$V3 == target.ind] <- IBS.ind$V2[IBS.ind$V3 == target.ind]
		IBS.rangelist[[i]] <- IBS.ind
	}
	print("rangelist completed")
	rangelists.chr.list[[j]] <- IBS.rangelist
	IBS.rangelist.h1 <- lapply(IBS.rangelist, function(x){x[x$target.hap == 1]})
	IBS.rangelist.h2 <- lapply(IBS.rangelist, function(x){x[x$target.hap == 2]})
	for(k in 1:length(bait.ns)){
		bait.samp <- bait[1:bait.ns[k]]
		for(l in 1:length(cM.cuts)){
			#print(l)
			cM.cut <- cM.cuts[l]			
			h1 <- reduce.in.bait(IBS.rangelist.h1, bait.samp, cM.cut)
			h2 <- reduce.in.bait(IBS.rangelist.h2, bait.samp, cM.cut)
			IBD2 <- mapply(FUN = function(g1,g2){ov <- findOverlaps(g1, g2) 
pintersect(g1[queryHits(ov)], g2[subjectHits(ov)])}, h1, h2)
			h1.covered[, j, k, l] <- sapply(h1, function(x){sum(width(x))})
			h2.covered[, j, k, l] <- sapply(h2, function(x){sum(width(x))})
			IBD2.covered[, j, k, l] <- sapply(IBD2, function(x){sum(width(x))})
			print(paste("chr", as.character(j), "bait", as.character(bait.ns[k]), "cM", as.character(cM.cut)))
			print(summary((h1.covered[, j, k, l] + h2.covered[, j, k, l])/2))
		}
	}

}

#total segments covered by any IBD on each chromosome
chr.lengths.ibd <- lapply(all.ranges.chr.list, function(x){sum(width(reduce(x)))})
tot.length <- sum(as.numeric(chr.lengths.ibd))

#save.image("figures/ibd_rangelists_872_lod3.RData")
#load("figures/ibd_rangelists_872_lod3.RData")

print(summary((h1.covered[, 1, 9, 1] + h2.covered[, 1, 9, 1])/2))



summary(apply(h1.covered[,,1,1], 1, sum))

avg.covered <- (h1.covered + h2.covered)/2

allchr.avg.covered <- apply(avg.covered, c(1,3,4), sum)
apply(allchr.avg.covered, c(2,3), median)


allchr.h1.covered <- apply(h1.covered, c(1,3,4), sum)
allchr.h2.covered <- apply(h2.covered, c(1,3,4), sum)
allchr.IBD2.covered <- apply(IBD2.covered, c(1,3,4), sum)

allchr.IBD1.covered <- allchr.avg.covered*2 - 2*allchr.IBD2.covered

pal <- c('#e66101','#fdb863','#b2abd2','#5e3c99')

#Given a matrix of medians with bait sizes in rows and cM cutoffs in columns,
#make a plot with bait sizes along horizontal axis
medians.plot <- function(median.mat, bait.ns, cM.cuts, pal, leg = FALSE, ...){
	plot(bait.ns, median.mat[,1], col = pal[1], type = "b", xlab = "Comparison sample size", bty = "n", las = 1, pch = 20, ...)
	for(i in 2:ncol(median.mat)){
		lines(bait.ns, median.mat[,i], col = pal[i], type = "b", pch = 20)
	}
	if(leg == TRUE){
		legend("topleft", legend = cM.cuts, pch = 20, col = pal, bty = "n", title = "cM cutoff")
	}
}

tiff("figures/median_cover_by_n_and_cM_lod3.tif", width = 9, height = 6, units = "in", res = 600, compression = "lzw")
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, .8, 0))
medians.plot(apply(allchr.avg.covered, c(2,3), median), bait.ns, cM.cuts, pal, leg = TRUE, ylab = "Median avg coverage (Gbp)", ylim = c(0, 1.6e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9), labels = c("0", ".5", "1", "1.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD2.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD2 coverage (Gbp)", ylim = c(0, 1.6e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9), labels = c("0", ".5", "1", "1.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD1.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD1 coverage (Gbp)", ylim = c(0, 1.6e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9), labels = c("0", ".5", "1", "1.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD1.covered + allchr.IBD2.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD1+ coverage (Gbp)", ylim = c(0, 1.6e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9), labels = c("0", ".5", "1", "1.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
dev.off()

pdf("figures/median_cover_by_n_and_cM_lod3.pdf", width = 9, height = 6)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.7, .7, 0), cex.lab = 1.2, cex.axis = 1.2)
medians.plot(apply(allchr.avg.covered, c(2,3), median), bait.ns, cM.cuts, pal, leg = TRUE, ylab = "Median avg coverage (Gbp)", ylim = c(0, 2.8e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9, 2e9, 2.5e9), labels = c("0", ".5", "1", "1.5", "2", "2.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD2.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD2 coverage (Gbp)", ylim = c(0, 2.8e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9, 2e9, 2.5e9), labels = c("0", ".5", "1", "1.5", "2", "2.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD1.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD1 coverage (Gbp)", ylim = c(0, 2.8e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9, 2e9, 2.5e9), labels = c("0", ".5", "1", "1.5", "2", "2.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
medians.plot(apply(allchr.IBD1.covered + allchr.IBD2.covered, c(2,3), median), bait.ns, cM.cuts, pal, ylab = "Median IBD1+ coverage (Gbp)", ylim = c(0, 2.8e9), yaxt = "n", xaxt = "n")
axis(2, at = c(0, 5e8, 1e9, 1.5e9, 2e9, 2.5e9), labels = c("0", ".5", "1", "1.5", "2", "2.5"), las = 1)
axis(1, at = c(100, 300, 500, 700, 872), las = 1)
dev.off()



#Numerical summaries of coverage reported in main text.

tot.length

apply(allchr.avg.covered, c(2,3), median)/tot.length
apply(allchr.IBD1.covered, c(2,3), median)/tot.length
apply(allchr.IBD2.covered, c(2,3), median)/tot.length
apply(allchr.IBD1.covered + allchr.IBD2.covered, c(2,3), median)/tot.length

apply(allchr.avg.covered, c(2,3), quantile, .9)/tot.length
apply(allchr.IBD1.covered, c(2,3), quantile, .9)/tot.length
apply(allchr.IBD2.covered, c(2,3), quantile, .9)/tot.length
apply(allchr.IBD1.covered + allchr.IBD2.covered, c(2,3), quantile, .9)/tot.length





