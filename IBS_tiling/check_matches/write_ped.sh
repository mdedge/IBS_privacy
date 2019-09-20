#!/bin/bash  

#9/20/19 by Doc Edge
#Goal: Check all the IBS calls on a chromosome for numbers of mismatches.
#Will write to a ped format and work with that.

chromo="19" #the chromosome to pull

#write to ped
vcftools --gzvcf ../phased/phased_HO_Euro_1kg_chr$chromo.vcf.gz --plink --chr $chromo --out phased_ped
#plink --vcf ../phased/phased_HO_Euro_1kg_chr$chromo.vcf.gz --recode --out phased_ped

#copy ibd calls into this folder
cp ../ibd/IBD_HO_Euro_1kg_chr$chromo.ibd ibd.ibd

#We do the rest in R
R

#Read in ped
ped <- read.table("phased_ped.ped", colClasses = "character")
#Read in map
map <- read.table("phased_ped.map", as.is = TRUE)
#Read in IBD calls
ibd <- read.table("ibd.ibd", as.is = TRUE)
ibd <- ibd[ibd$V9 > 1,]

#Reformat genotypes of ped table so that there is one haplotype in each row
hap.id <- paste(rep(ped$V1, each = 2), rep(c("1", "2"), length(ped$V1)), sep = ".")

haps.1 <- as.matrix(ped[,seq(7, ncol(ped) - 1, 2)])
haps.2 <- as.matrix(ped[,seq(8, ncol(ped), 2)])

haps <- rbind(haps.1, haps.2)[order(sequence( c(nrow(haps.1), nrow(haps.2) )  )), ]

ibd$hapA <- paste(ibd$V1, ibd$V2, sep = ".")
ibd$hapB <- paste(ibd$V3, ibd$V4, sep = ".")

#Compute matches for each shared segment

nsnps <- numeric(nrow(ibd))
nmatch <- numeric(nrow(ibd))
nmatch.altA <- numeric(nrow(ibd))
nmatch.altB <- numeric(nrow(ibd))
nmatch.altAB <- numeric(nrow(ibd))
nmatch.max <- numeric(nrow(ibd))
nmatch.randoA <- numeric(nrow(ibd))
nmatch.rando2 <- numeric(nrow(ibd))

for(i in 1:nrow(ibd)){
	#extract the two segments
	start <- ibd$V6[i]
	end <- ibd$V7[i]
	cols <- which(map$V4 > start & map$V4 < end)
	haploA <- haps[hap.id == ibd$hapA[i], cols]
	haploA.alt <- haps[hap.id == paste(ibd$V1[i], 3 - ibd$V2[i], sep = "."), cols]
	haploB <- haps[hap.id == ibd$hapB[i], cols]
	haploB.alt <- haps[hap.id == paste(ibd$V3[i], 3 - ibd$V4[i], sep = "."), cols]
	randohap1 <- haps[sample(1:nrow(haps),1), cols]
	randohap2 <- haps[sample(1:nrow(haps),1), cols]
	nsnps[i] <- length(cols)
	nmatch[i] <- sum(haploA == haploB)
	nmatch.altA[i] <- sum(haploA.alt == haploB)
	nmatch.altB[i] <- sum(haploA == haploB.alt)
	nmatch.altAB[i] <- sum(haploA.alt == haploB.alt)
	nmatch.max[i] <- max(c(nmatch[i], nmatch.altA[i], nmatch.altB[i], nmatch.altAB[i]))
	nmatch.randoA[i] <- sum(haploA == randohap1)
	nmatch.rando2[i] <- sum(randohap1 == randohap2)
	if(i %% 10000 == 0){print(i)}
}

summary(nmatch / nsnps)
summary(nmatch.altA / nsnps)
summary(nmatch.altB / nsnps)
summary(nmatch.altAB / nsnps)
summary(nmatch.max / nsnps)
summary(nmatch.randoA / nsnps)
summary(nmatch.rando2 / nsnps)

#plot lod score vs. haplotype accuracy
plot(ibd$V8, nmatch / nsnps, xlab = "LOD score", xlim = c(0, 20))

plot(ibd$V8, nmatch.max / nsnps, xlab = "LOD score", xlim = c(0, 20))


jpeg("LOD_and_matchrate.jpg")
plot(ibd$V8, nmatch.max / nsnps, xlab = "LOD score", xlim = c(0, 20), pch = ".", ylab = "match rate", bty = "n")
dev.off()

summary(nmatch.max / nsnps)
summary((nmatch.max / nsnps)[ibd$V8 > 2])
summary((nmatch.max / nsnps)[ibd$V8 > 3])
summary((nmatch.max / nsnps)[ibd$V8 > 4])
summary((nmatch.max / nsnps)[ibd$V8 > 5])
summary((nmatch.max / nsnps)[ibd$V8 > 10])



#Compute compatibility for shared segments

n.incompat <- numeric(nrow(ibd))

person.hap <- rep(ped$V1, each = 2)
#function to identify incompatible homozygotes
incompatible.hom <- function(x){
	if(x[1] == x[2] & x[3] == x[4] & x[1] != x[3]){
		return(TRUE)
	}
	FALSE
}

for(i in 1:nrow(ibd)){
	#extract the two segments
	start <- ibd$V6[i]
	end <- ibd$V7[i]
	cols <- which(map$V4 > start & map$V4 < end)
	haplosA <- haps[person.hap == ibd$V1[i], cols]
	haplosB <- haps[person.hap == ibd$V3[i], cols]

	n.incompat[i] <- sum(apply(rbind(haplosA, haplosB), 2, incompatible.hom ))
	if(i %% 10000 == 0){print(i)}
}



