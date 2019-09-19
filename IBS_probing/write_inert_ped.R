

#Goal: sample random alleles and write them in ped format

dat <- read.table("allelefreq.frq", head = FALSE, skip = 1, stringsAsFactors = FALSE)

#extract alleles and frequencies
A1 <- sapply(strsplit(dat$V5, ":"), "[[", 1)
A1.frq <- as.numeric(sapply(strsplit(dat$V5, ":"), "[[", 2))
A2 <- sapply(strsplit(dat$V6, ":"), "[[", 1)


#To prevent IBS matches in low-diversity regions, artificially inflate rare allele
#frequencies to .1
A1.frq[A1.frq < .1] <- .1
A1.frq[A1.frq > .9] <- .9

#sample alleles according to their frequencies
set.seed(8675309)
saf.pick.min <- rbinom(nrow(dat) * max(dat$V4), 1, rep(1-A1.frq, max(dat$V4) ) )
samp.maf <- character(nrow(dat) * max(dat$V4))

A1.rep <- rep(A1, max(dat$V4))
A2.rep <- rep(A2, max(dat$V4))
samp.maf[saf.pick.min == 0] <- as.character(A2.rep[saf.pick.min == 0])
samp.maf[saf.pick.min == 1] <- as.character(A1.rep[saf.pick.min == 1])

#format for ped
mat1 <- matrix(samp.maf[1:(length(samp.maf)/2)], byrow = TRUE, ncol = nrow(dat))
mat2 <- matrix(samp.maf[(length(samp.maf)/2 + 1):(length(samp.maf))], byrow = TRUE, ncol = nrow(dat))

l <- list(mat1, mat2)
geno.mat <- do.call(cbind, l)[, order(sequence(sapply(l, ncol)))]

#make a ped file
ped <- read.table("phased_ped.ped", head = FALSE, stringsAsFactors = FALSE)

ped.inert <- cbind(ped[,1:6], geno.mat)

write.table( ped.inert, "IBSinert.ped", quote = FALSE, row.names = FALSE, col.names = FALSE)


