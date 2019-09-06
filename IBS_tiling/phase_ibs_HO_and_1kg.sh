#!/bin/bash

#9/4/19 by Doc Edge
#Goal: Phase and call IBS to reproduce the IBS tiling results of Edge and Coop, "Attacks on genetic 
#privacy via uploads to genealogical databases." 
#This script requires 

#### java (to run Beagle and refinedIBD)
#### beagle version 5 available (provided)
#### refinedIBD available (provided)


#########################################################



#Beagle 5 phasing
date
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=22 out=phased/phased_HO_Euro_1kg_chr22 map=plink.GRCh37.map/plink.chr22.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=21 out=phased/phased_HO_Euro_1kg_chr21 map=plink.GRCh37.map/plink.chr21.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=20 out=phased/phased_HO_Euro_1kg_chr20 map=plink.GRCh37.map/plink.chr20.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=19 out=phased/phased_HO_Euro_1kg_chr19 map=plink.GRCh37.map/plink.chr19.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=18 out=phased/phased_HO_Euro_1kg_chr18 map=plink.GRCh37.map/plink.chr18.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=17 out=phased/phased_HO_Euro_1kg_chr17 map=plink.GRCh37.map/plink.chr17.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=16 out=phased/phased_HO_Euro_1kg_chr16 map=plink.GRCh37.map/plink.chr16.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=15 out=phased/phased_HO_Euro_1kg_chr15 map=plink.GRCh37.map/plink.chr15.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=14 out=phased/phased_HO_Euro_1kg_chr14 map=plink.GRCh37.map/plink.chr14.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=13 out=phased/phased_HO_Euro_1kg_chr13 map=plink.GRCh37.map/plink.chr13.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=12 out=phased/phased_HO_Euro_1kg_chr12 map=plink.GRCh37.map/plink.chr12.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=11 out=phased/phased_HO_Euro_1kg_chr11 map=plink.GRCh37.map/plink.chr11.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=10 out=phased/phased_HO_Euro_1kg_chr10 map=plink.GRCh37.map/plink.chr10.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=9 out=phased/phased_HO_Euro_1kg_chr9 map=plink.GRCh37.map/plink.chr9.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=8 out=phased/phased_HO_Euro_1kg_chr8 map=plink.GRCh37.map/plink.chr8.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=7 out=phased/phased_HO_Euro_1kg_chr7 map=plink.GRCh37.map/plink.chr7.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=6 out=phased/phased_HO_Euro_1kg_chr6 map=plink.GRCh37.map/plink.chr6.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=5 out=phased/phased_HO_Euro_1kg_chr5 map=plink.GRCh37.map/plink.chr5.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=4 out=phased/phased_HO_Euro_1kg_chr4 map=plink.GRCh37.map/plink.chr4.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=3 out=phased/phased_HO_Euro_1kg_chr3 map=plink.GRCh37.map/plink.chr3.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=2 out=phased/phased_HO_Euro_1kg_chr2 map=plink.GRCh37.map/plink.chr2.GRCh37.map nthreads=4 
java -jar beagle.16May19.351.jar gt=HO_Euro_1kg_biallelic.vcf.gz chrom=1 out=phased/phased_HO_Euro_1kg_chr1 map=plink.GRCh37.map/plink.chr1.GRCh37.map nthreads=4 
date


#Use refinedIBD to call IBD on phased vcf
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr1.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr1 map=plink.GRCh37.map/plink.chr1.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr2.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr2 map=plink.GRCh37.map/plink.chr2.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr3.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr3 map=plink.GRCh37.map/plink.chr3.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr4.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr4 map=plink.GRCh37.map/plink.chr4.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr5.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr5 map=plink.GRCh37.map/plink.chr5.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr6.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr6 map=plink.GRCh37.map/plink.chr6.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr7.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr7 map=plink.GRCh37.map/plink.chr7.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr8.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr8 map=plink.GRCh37.map/plink.chr8.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr9.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr9 map=plink.GRCh37.map/plink.chr9.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr10.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr10 map=plink.GRCh37.map/plink.chr10.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr11.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr11 map=plink.GRCh37.map/plink.chr11.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr12.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr12 map=plink.GRCh37.map/plink.chr12.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr13.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr13 map=plink.GRCh37.map/plink.chr13.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr14.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr14 map=plink.GRCh37.map/plink.chr14.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr15.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr15 map=plink.GRCh37.map/plink.chr15.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr16.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr16 map=plink.GRCh37.map/plink.chr16.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr17.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr17 map=plink.GRCh37.map/plink.chr17.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr18.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr18 map=plink.GRCh37.map/plink.chr18.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr19.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr19 map=plink.GRCh37.map/plink.chr19.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr20.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr20 map=plink.GRCh37.map/plink.chr20.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr21.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr21 map=plink.GRCh37.map/plink.chr21.GRCh37.map nthreads=4 length=.8 lod=1
java -Xss5m -Xmx50g -jar refined-ibd.16May19.ad5.jar gt=phased/phased_HO_Euro_1kg_chr22.vcf.gz out=ibd/IBD_HO_Euro_1kg_chr22 map=plink.GRCh37.map/plink.chr22.GRCh37.map nthreads=4 length=.8 lod=1

gunzip ibd/IBD_HO_Euro_1kg_chr1.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr2.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr3.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr4.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr5.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr6.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr7.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr8.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr9.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr10.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr11.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr12.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr13.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr14.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr15.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr16.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr17.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr18.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr19.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr20.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr21.ibd
gunzip ibd/IBD_HO_Euro_1kg_chr22.ibd




