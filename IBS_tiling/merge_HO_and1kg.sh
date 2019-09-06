

#9/4/19 by Doc Edge
#Goal: Prepare data to reproduce the IBS tiling results of Edge and Coop, "Attacks on genetic 
#privacy via uploads to genealogical databases." 
#The goal of this script is only to produce the file HO_Euro_1kg.vcf.gz from the Human Origins
#data and the 1000 Genomes data. It can be skipped if the user would rather just use
#the provided HO_Euro_1kg.vcf.gz data file; it is provided mainly so that we can
#show our steps in creating this dataset.

#This script requires
#### plink version 1.9
#### R 
#### vcftools
#### bcftools
#### eigensoft (to run convertf to convert Human Origins data from eigenstrat to ped)

# We also assume that the human origins and 1000 Genomes phase 3 data are 
# available on the local machine. 
#In this script, the directory in which Human Origins data sit is /home/shared/datasets/human_origins/
#The 1000 Genomes directory is /home/medge/Data/1000gen_rel_20130502/
#You will need to replace these directory names to match your system in this file and in the provided
#par.PED.EIGENSTRAT file.

#The Human Origins data were downloaded from https://reich.hms.harvard.edu/downloadable-genotypes-present-day-and-ancient-dna-data-compiled-published-papers . We used the 1240K+HO dataset, version 37.2

#The 1000 Genomes data were downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/


##############################################################


#Convert Human Origins (eigenstrat format) to ped so that we can get site positions easily
#(Assumes Eigensoft is a directory w/in user's home directory)
#~/EIG-7.2.1/src/convertf -p par.PED.EIGENSTRAT


#pull the European sample ids from 1000 genomes
grep EUR /home/medge/Data/1000gen_rel_20130502/integrated_call_samples_v3.20130502.ALL.panel | cut -f1 > out/EUR.sample.ids


#########################################################

#Select individuals from Human Origins. 
#In particular, want modern Europeans
#and want to pull in 1000Genomes data ourselves to avoid using the haploid 
#".SG" entries in HO.

R

dat <- read.delim("/home/shared/datasets/human_origins/v37.2.1240K_HumanOrigins.clean4.txt",head=TRUE,sep="\t",as.is=TRUE)

pops <- unique(dat[dat[,6] < 1,4])
pops[!grepl(".DG", pops) & !grepl(".SG", pops) & !grepl("Ignore", pops)]
eur.grps <- c("French", "Sardinian", "Orcadian", "Russian", "Italian_North", "Italian_South", "Basque", "Bulgarian", "Hungarian","Lithuanian", "Ukrainian", "Estonian", "Belarusian", "Mordovian", "Czech", "Icelandic", "Greek", "Scottish", "English", "Spanish", "Spanish_North", "Finnish", "Maltese", "Croatian", "Norwegian", "Sicilian", "Albanian",  "Romanian", "Adygei")

summary(as.factor(dat[dat[,4] %in% eur.grps,4]))

incl.grp <-dat[,4] %in% c(eur.grps)
dim(dat[incl.grp, c(2,4)])

#Select only the modern Europeans who are not in 1000 Genomes 
incl.grp.no1kg <- incl.grp - grepl("HG0", dat[,2])
#incl.grp.no1kg <- incl.grp - (dat[,5] == "1KGPhase3") #This doesn't work. 74 1kg samples not marked as coming from 1kg.
incl.grp.no1kg[incl.grp.no1kg == -1] <- 0
incl.grp.no1kg <- as.logical(incl.grp.no1kg)

summary(as.factor(dat[incl.grp.no1kg & as.numeric(dat[,6]) < 1,4]))

dat2 <- read.table("/home/shared/datasets/human_origins/v37.2.1240K_HumanOrigins.pedind")

ids.incl <- dat2[dat2$V2 %in% dat[incl.grp.no1kg,2],1:2]

#Check that filtering out 1kg was successful
ids.1kg <- read.table("out/EUR.sample.ids", stringsAsFactors = FALSE)$V1
sum(ids.incl[,2] %in% ids.1kg)

write.table(ids.incl, "out/Euro_inclgrps_inds_HO_plink.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

q()
n

######################################################



#Get the site positions for Human Origins so we can pull those from 1000 Genomes (both on genome build 37)
awk -v OFS='\t' '{print $1, $4}' /home/shared/datasets/human_origins/v37.2.1240K_HumanOrigins.map > out/HO_sitepos.txt



#Next, select just EUR samples from 1000 Genomes and filter down to Human Origins sites.
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr1
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr2
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr3
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr4
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr5
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr6
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr7
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr8
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr9
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr10
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr11
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr12
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr13
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr14
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr15
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr16
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr17
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr18
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr19
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr20
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr21
vcftools --gzvcf /home/medge/Data/1000gen_rel_20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep out/EUR.sample.ids --positions out/HO_sitepos.txt --recode --out out/EUR_1kg_HOpos_chr22

#Concatenate chromosomes into one file
echo "out/EUR_1kg_HOpos_chr1.recode.vcf
out/EUR_1kg_HOpos_chr2.recode.vcf
out/EUR_1kg_HOpos_chr3.recode.vcf
out/EUR_1kg_HOpos_chr4.recode.vcf
out/EUR_1kg_HOpos_chr5.recode.vcf
out/EUR_1kg_HOpos_chr6.recode.vcf
out/EUR_1kg_HOpos_chr7.recode.vcf
out/EUR_1kg_HOpos_chr8.recode.vcf
out/EUR_1kg_HOpos_chr9.recode.vcf
out/EUR_1kg_HOpos_chr10.recode.vcf
out/EUR_1kg_HOpos_chr11.recode.vcf
out/EUR_1kg_HOpos_chr12.recode.vcf
out/EUR_1kg_HOpos_chr13.recode.vcf
out/EUR_1kg_HOpos_chr14.recode.vcf
out/EUR_1kg_HOpos_chr15.recode.vcf
out/EUR_1kg_HOpos_chr16.recode.vcf
out/EUR_1kg_HOpos_chr17.recode.vcf
out/EUR_1kg_HOpos_chr18.recode.vcf
out/EUR_1kg_HOpos_chr19.recode.vcf
out/EUR_1kg_HOpos_chr20.recode.vcf
out/EUR_1kg_HOpos_chr21.recode.vcf
out/EUR_1kg_HOpos_chr22.recode.vcf" >> out/merge_fns.txt

bcftools concat --file-list out/merge_fns.txt --output-type z --output out/EUR_1kg_HOpos_allchr.gzvcf

#Remove non-SNPs
vcftools --gzvcf out/EUR_1kg_HOpos_allchr.gzvcf --remove-indels --recode --recode-INFO-all --out out/EUR_1kg_HOpos_allchr_SNPs_only

#Make a file recording the reference / alt alleles for each locus included.
awk '{ print $1,  $2,  $3,  $4,  $5 }' out/EUR_1kg_HOpos_allchr_SNPs_only.recode.vcf > out/ref_rsid_1kg_HOpos.txt
sed '/^#/ d' < out/ref_rsid_1kg_HOpos.txt > out/ref_rsid_1kg_HOpos_nohead.txt
awk '{print $3}' < out/ref_rsid_1kg_HOpos_nohead.txt > out/rsid_1kg.txt


########################################################

#Convert Human Origins  to bed with AGCT alleles, and then to vcf.
plink --file /home/shared/datasets/human_origins/v37.2.1240K_HumanOrigins --make-bed --alleleACGT --keep out/Euro_inclgrps_inds_HO_plink.txt --out out/HO_Euro_AGCT --allow-no-sex --extract out/rsid_1kg.txt

#Assign the ref alleles to match 1000Genomes
plink --bfile out/HO_Euro_AGCT --recode vcf-iid --a2-allele out/ref_rsid_1kg_HOpos_nohead.txt 4 3 --out out/HO_Euro_vcf --allow-no-sex


#Now filter the 1kg again to contain only the SNPs that also appear in HO.
awk '{print $1,  $2,  $3}' < out/HO_Euro_vcf.vcf > out/pos_1kg_and_HO.txt
sed '/^#/ d' < out/pos_1kg_and_HO.txt > out/pos_1kg_and_HO_nohead.txt
awk '{print $3}' < out/pos_1kg_and_HO_nohead.txt > out/rsid_1kg_HO.txt

vcftools --vcf out/EUR_1kg_HOpos_allchr_SNPs_only.recode.vcf --snps out/rsid_1kg_HO.txt --recode --out out/EUR_1kg_HOpos_allchr_inHO.vcf


#Zip to prepare for merge
bgzip out/HO_Euro_vcf.vcf
bgzip out/EUR_1kg_HOpos_allchr_inHO.vcf.recode.vcf

#Human Origins and 1kg occasionally use different SNP IDs, so we write new IDs so that they match
#As of 9/4/19, this step doesn't matter bc we've filtered to only sites with matching IDs
bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' out/EUR_1kg_HOpos_allchr_inHO.vcf.recode.vcf.gz -o out/EUR_1kg_HOpos_allchr_SNPs_only.recode.IDcpra.vcf.gz
bcftools annotate -Oz -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' out/HO_Euro_vcf.vcf.gz  -o out/HO_Euro_vcf.IDcpra.vcf.gz

#make indices to prepare for merge
bcftools index out/HO_Euro_vcf.IDcpra.vcf.gz
bcftools index out/EUR_1kg_HOpos_allchr_SNPs_only.recode.IDcpra.vcf.gz

vcftools --gzvcf out/HO_Euro_vcf.IDcpra.vcf.gz
vcftools --gzvcf out/EUR_1kg_HOpos_allchr_SNPs_only.recode.IDcpra.vcf.gz

#merge HO and 1000 genomes
bcftools merge -m both -O z -o out/HO_Euro_1kg.vcf.gz out/HO_Euro_vcf.IDcpra.vcf.gz out/EUR_1kg_HOpos_allchr_SNPs_only.recode.IDcpra.vcf.gz

#vcftools --gzvcf out/HO_Euro_1kg.vcf.gz 
bcftools view -m2 -M2 -v snps out/HO_Euro_1kg.vcf.gz -Oz -o out/HO_Euro_1kg_biallelic.vcf.gz

#vcftools --gzvcf out/HO_Euro_1kg_biallelic.vcf.gz

#Copy the merged file into the main directory. 
cp out/HO_Euro_1kg_biallelic.vcf.gz HO_Euro_1kg_biallelic.vcf.gz 






