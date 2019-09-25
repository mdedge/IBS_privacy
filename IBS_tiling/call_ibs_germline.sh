#!/bin/bash

#9/22/19 by Doc Edge
#Goal: Use germline to call IBS for Edge and Coop, "Attacks on genetic 
#privacy via uploads to genealogical databases." 
#This script requires 

#### germline


#########################################################

#~/germline-1-5-3/germline -min_m 1cM -haploid -input /home/shared/datasets/human_origins/phase_070119/phased_HO_EuroC_1kg_chr1.ped genetic_maps/phased_HO_EuroC_1kg_complete_1.map  -output germline_ibd/HO_EuroC_1kg_chr_1_cM_haploid &
#My go at baiting 
#~/germline-1-5-3/germline -min_m 1cM -bits 64 -w_extend -g_extend -input /home/shared/datasets/human_origins/phase_070119/bait_hets_phased_HO_EuroC_1kg_chr19.ped genetic_maps/phased_HO_EuroC_1kg_complete_19.map  -output germline_ibd/het_bait &



#First, convert the phased vcfs to phased ped.

vcftools --gzvcf phased/phased_HO_Euro_1kg_chr1.vcf.gz --plink --chr 1 --out phased_ped/HO_Euro_1kg_chr1
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr2.vcf.gz --plink --chr 2 --out phased_ped/HO_Euro_1kg_chr2
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr3.vcf.gz --plink --chr 3 --out phased_ped/HO_Euro_1kg_chr3
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr4.vcf.gz --plink --chr 4 --out phased_ped/HO_Euro_1kg_chr4
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr5.vcf.gz --plink --chr 5 --out phased_ped/HO_Euro_1kg_chr5
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr6.vcf.gz --plink --chr 6 --out phased_ped/HO_Euro_1kg_chr6
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr7.vcf.gz --plink --chr 7 --out phased_ped/HO_Euro_1kg_chr7
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr8.vcf.gz --plink --chr 8 --out phased_ped/HO_Euro_1kg_chr8
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr9.vcf.gz --plink --chr 9 --out phased_ped/HO_Euro_1kg_chr9
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr10.vcf.gz --plink --chr 10 --out phased_ped/HO_Euro_1kg_chr10
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr11.vcf.gz --plink --chr 11 --out phased_ped/HO_Euro_1kg_chr11
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr12.vcf.gz --plink --chr 12 --out phased_ped/HO_Euro_1kg_chr12
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr13.vcf.gz --plink --chr 13 --out phased_ped/HO_Euro_1kg_chr13
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr14.vcf.gz --plink --chr 14 --out phased_ped/HO_Euro_1kg_chr14
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr15.vcf.gz --plink --chr 15 --out phased_ped/HO_Euro_1kg_chr15
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr16.vcf.gz --plink --chr 16 --out phased_ped/HO_Euro_1kg_chr16
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr17.vcf.gz --plink --chr 17 --out phased_ped/HO_Euro_1kg_chr17
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr18.vcf.gz --plink --chr 18 --out phased_ped/HO_Euro_1kg_chr18
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr19.vcf.gz --plink --chr 19 --out phased_ped/HO_Euro_1kg_chr19
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr20.vcf.gz --plink --chr 20 --out phased_ped/HO_Euro_1kg_chr20
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr21.vcf.gz --plink --chr 21 --out phased_ped/HO_Euro_1kg_chr21
vcftools --gzvcf phased/phased_HO_Euro_1kg_chr22.vcf.gz --plink --chr 22 --out phased_ped/HO_Euro_1kg_chr22


#make copies of the chromosome maps that linearly interpolate all the sites we have
R

for(i in 1:22){
	cm.map <- read.table(paste("plink.GRCh37.map/plink.chr", i, ".GRCh37.map", sep = ""))
	site.map <- read.table(paste("phased_ped/HO_Euro_1kg_chr", i, ".map", sep = ""))
	interp.cm <- approx(cm.map$V4, cm.map$V3, xout = site.map$V4, rule = 2)$y
	site.map$V3 <- interp.cm 
	write.table(site.map, paste("plink.GRCh37.map.allsites/plink.chr", i, ".GRCh37.allsites.map", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE )
}



q()
n


#Call ibd with haploid flag
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr1.ped plink.GRCh37.map.allsites/plink.chr1.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr1_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr2.ped plink.GRCh37.map.allsites/plink.chr2.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr2_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr3.ped plink.GRCh37.map.allsites/plink.chr3.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr3_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr4.ped plink.GRCh37.map.allsites/plink.chr4.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr4_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr5.ped plink.GRCh37.map.allsites/plink.chr5.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr5_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr6.ped plink.GRCh37.map.allsites/plink.chr6.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr6_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr7.ped plink.GRCh37.map.allsites/plink.chr7.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr7_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr8.ped plink.GRCh37.map.allsites/plink.chr8.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr8_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr9.ped plink.GRCh37.map.allsites/plink.chr9.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr9_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr10.ped plink.GRCh37.map.allsites/plink.chr10.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr10_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr11.ped plink.GRCh37.map.allsites/plink.chr11.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr11_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr12.ped plink.GRCh37.map.allsites/plink.chr12.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr12_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr13.ped plink.GRCh37.map.allsites/plink.chr13.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr13_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr14.ped plink.GRCh37.map.allsites/plink.chr14.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr14_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr15.ped plink.GRCh37.map.allsites/plink.chr15.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr15_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr16.ped plink.GRCh37.map.allsites/plink.chr16.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr16_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr17.ped plink.GRCh37.map.allsites/plink.chr17.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr17_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr18.ped plink.GRCh37.map.allsites/plink.chr18.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr18_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr19.ped plink.GRCh37.map.allsites/plink.chr19.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr19_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr20.ped plink.GRCh37.map.allsites/plink.chr20.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr20_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr21.ped plink.GRCh37.map.allsites/plink.chr21.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr21_haploid -err_hom 0 -err_het 0 -w_extend -bits 32
./germline -min_m 1cM -haploid -input phased_ped/HO_Euro_1kg_chr22.ped plink.GRCh37.map.allsites/plink.chr22.GRCh37.allsites.map -output ibd_germline_haploid/ibd_chr22_haploid -err_hom 0 -err_het 0 -w_extend -bits 32



#Call ibd with h_extend flag
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr1.ped plink.GRCh37.map.allsites/plink.chr1.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr1_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr2.ped plink.GRCh37.map.allsites/plink.chr2.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr2_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr3.ped plink.GRCh37.map.allsites/plink.chr3.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr3_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr4.ped plink.GRCh37.map.allsites/plink.chr4.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr4_hextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr5.ped plink.GRCh37.map.allsites/plink.chr5.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr5_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr6.ped plink.GRCh37.map.allsites/plink.chr6.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr6_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr7.ped plink.GRCh37.map.allsites/plink.chr7.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr7_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr8.ped plink.GRCh37.map.allsites/plink.chr8.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr8_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr9.ped plink.GRCh37.map.allsites/plink.chr9.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr9_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr10.ped plink.GRCh37.map.allsites/plink.chr10.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr10_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr11.ped plink.GRCh37.map.allsites/plink.chr11.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr11_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr12.ped plink.GRCh37.map.allsites/plink.chr12.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr12_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr13.ped plink.GRCh37.map.allsites/plink.chr13.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr13_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr14.ped plink.GRCh37.map.allsites/plink.chr14.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr14_hextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr15.ped plink.GRCh37.map.allsites/plink.chr15.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr15_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr16.ped plink.GRCh37.map.allsites/plink.chr16.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr16_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr17.ped plink.GRCh37.map.allsites/plink.chr17.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr17_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr18.ped plink.GRCh37.map.allsites/plink.chr18.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr18_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr19.ped plink.GRCh37.map.allsites/plink.chr19.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr19_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr20.ped plink.GRCh37.map.allsites/plink.chr20.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr20_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr21.ped plink.GRCh37.map.allsites/plink.chr21.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr21_hextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -h_extend -input phased_ped/HO_Euro_1kg_chr22.ped plink.GRCh37.map.allsites/plink.chr22.GRCh37.allsites.map -output ibd_germline_hextend/ibd_chr22_hextend -err_hom 0 -err_het 0 -w_extend 



#Call ibd with g_extend flag (i.e. looking for incompatible homozygotes)
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr1.ped plink.GRCh37.map.allsites/plink.chr1.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr1_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr2.ped plink.GRCh37.map.allsites/plink.chr2.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr2_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr3.ped plink.GRCh37.map.allsites/plink.chr3.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr3_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr4.ped plink.GRCh37.map.allsites/plink.chr4.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr4_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr5.ped plink.GRCh37.map.allsites/plink.chr5.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr5_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr6.ped plink.GRCh37.map.allsites/plink.chr6.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr6_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr7.ped plink.GRCh37.map.allsites/plink.chr7.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr7_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr8.ped plink.GRCh37.map.allsites/plink.chr8.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr8_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr9.ped plink.GRCh37.map.allsites/plink.chr9.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr9_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr10.ped plink.GRCh37.map.allsites/plink.chr10.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr10_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr11.ped plink.GRCh37.map.allsites/plink.chr11.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr11_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr12.ped plink.GRCh37.map.allsites/plink.chr12.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr12_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr13.ped plink.GRCh37.map.allsites/plink.chr13.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr13_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr14.ped plink.GRCh37.map.allsites/plink.chr14.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr14_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr15.ped plink.GRCh37.map.allsites/plink.chr15.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr15_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr16.ped plink.GRCh37.map.allsites/plink.chr16.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr16_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr17.ped plink.GRCh37.map.allsites/plink.chr17.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr17_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr18.ped plink.GRCh37.map.allsites/plink.chr18.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr18_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr19.ped plink.GRCh37.map.allsites/plink.chr19.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr19_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr20.ped plink.GRCh37.map.allsites/plink.chr20.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr20_gextend -err_hom 0 -err_het 0 -w_extend 
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr21.ped plink.GRCh37.map.allsites/plink.chr21.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr21_gextend -err_hom 0 -err_het 0 -w_extend
./germline -min_m 1cM -g_extend -input phased_ped/HO_Euro_1kg_chr22.ped plink.GRCh37.map.allsites/plink.chr22.GRCh37.allsites.map -output ibd_germline_gextend/ibd_chr22_gextend -err_hom 0 -err_het 0 -w_extend




