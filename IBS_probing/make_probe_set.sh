#!/bin/bash  

#9/16/19 by Doc Edge
#Goal: Produce a dataset of IBS probes for a given chromosome and region by filling in
#genotypes outside the region with randomly drawn alleles, then merge it with the real
#data and call IBD. Requires plink, java, bcftools, and vcftools.

#Steps are:
#1) extract the chromosome and phase it.
#2) extract just the probe part
#3) compute allele frequencies
#4) Make an IBS-inert version of the file where alleles are filled in randomly at their frequencies
#5) cut the parts other than the probing region from the inert file and splice them together
#with the probe part from the real file.
#6) Merge the real data with the probe dataset and run refinedIBD on the whole thing.

#Note: It'd save a bit of time just to phase the probe part, but this doesn't take very long as it is.

chromo="19" #the chromosome to pull
lbcoord="45136023" #the lower bound of the probe part
#lbcoord="43834377" #the lower bound of the probe part
ubcoord="45817621" #the upper bound of the probe part
#ubcoord="47047801" #the upper bound of the probe part
lbcl=$lbcoord-1
ubcr=$ubcoord+1

#Extract and phase the dataset
java -jar ../IBS_tiling/beagle.16May19.351.jar gt=../IBS_tiling/HO_Euro_1kg_biallelic.vcf.gz chrom=$chromo out=phased_HO_Euro_1kg_chr$chromo map=../IBS_tiling/plink.GRCh37.map/plink.chr$chromo.GRCh37.map nthreads=4 
#write this to ped
plink --vcf phased_HO_Euro_1kg_chr$chromo.vcf.gz --recode --out phased_ped

#Use vcftools to pull out just the active part.
vcftools --gzvcf phased_HO_Euro_1kg_chr$chromo.vcf.gz --chr $chromo --from-bp $lbcoord --to-bp $ubcoord --recode --out EUR_1kg_chr$chromo.active

#Compute allele frequencies across the chromosome
vcftools --gzvcf phased_HO_Euro_1kg_chr$chromo.vcf.gz --freq --out allelefreq 


#Make a .ped-format file with sampled alleles at all sites
Rscript write_inert_ped.R
cp phased_ped.map IBSinert.map
cp phased_ped.nosex IBSinert.nosex

#Make the reference alleles match those in the original vcf so that we can merge files later

####first pull out the ref allele information
gunzip phased_HO_Euro_1kg_chr$chromo.vcf.gz
awk '{ print $1,  $2,  $3,  $4,  $5 }' phased_HO_Euro_1kg_chr$chromo.vcf > ref_rsid.txt
sed '/^#/ d' < ref_rsid.txt > ref_rsid_nohead.txt
bgzip phased_HO_Euro_1kg_chr$chromo.vcf

####then make plink match the reference alleles and convert to vcf
plink --file IBSinert --recode vcf-iid --a2-allele ref_rsid_nohead.txt 4 3 --out IBSinert --allow-no-sex


#pull out just the lower and upper inert parts
vcftools --vcf IBSinert.vcf --chr $chromo --to-bp $lbcl --recode --out left_inert
vcftools --vcf IBSinert.vcf --chr $chromo --from-bp $ubcr --recode --out right_inert

#Splice together
bcftools concat left_inert.recode.vcf EUR_1kg_chr$chromo.active.recode.vcf --output leftprobe.vcf --output-type v
bcftools concat leftprobe.vcf right_inert.recode.vcf --output probe.vcf --output-type v

#change individual ids to have a .probe attached.
bcftools query -l probe.vcf > sampleids.txt
sed -e 's/$/.probe/' "sampleids.txt" > sampleids.probe.txt
bcftools reheader -s sampleids.probe.txt probe.vcf --output probe.head.vcf 

#make the inert parts look as though they're phased so they can be fed to refinedIBD
#i.e. replace / with |
sed 's/\//|/g' probe.head.vcf > probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf

#Clean up
rm leftprobe.vcf
rm allelefreq.*
rm EUR_1kg_chr$chromo.*
rm IBSinert.*
rm left_inert.*
rm phased_ped.*
rm right_inert.*
rm sampleids.*
rm probe.vcf
rm probe.head.vcf
rm ref_rsid*

#Make indices to prep for merge
bgzip probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf
bcftools index probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz
bcftools index phased_HO_Euro_1kg_chr$chromo.vcf.gz

#vcftools --gzvcf probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz --freq --out allelefreq_probe 

bcftools merge -m both -O z -o real_and_probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz phased_HO_Euro_1kg_chr$chromo.vcf.gz probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz


#Run refined IBD on this file
java -Xss5m -Xmx50g -jar ../IBS_tiling/refined-ibd.16May19.ad5.jar gt=real_and_probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz out=IBD_real_and_probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz map=../IBS_tiling/plink.GRCh37.map/plink.chr$chromo.GRCh37.map nthreads=4 length=.8 lod=1

gunzip IBD_real_and_probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz.ibd.gz


#################3
#phase with germline in haploid mode

#First, convert file to ped
plink --vcf real_and_probe.chr$chromo.lb$lbcoord.ub$ubcoord.vcf.gz --recode --out probe_ped.lb$lbcoord.ub$ubcoord

#call IBD with germline
../IBS_tiling/germline -min_m 1cM -haploid -input probe_ped.lb$lbcoord.ub$ubcoord.ped ../IBS_tiling/plink.GRCh37.map.allsites/plink.chr$chromo.GRCh37.allsites.map -output ibd_rp_chr$chromo.lb$lbcoord.ub$ubcoord.germlinehaploid -err_hom 0 -err_het 0 -w_extend -bits 16

rm *.map
rm *.nosex
rm *.log


