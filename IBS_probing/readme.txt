Scripts to run the IBD probing analysis in Edge & Coop, "Attacks on Genetic Privacy via Uploads to Genealogical Databases."

To carry out the analysis, one needs to be able to run R, java, and germline from the terminal.

Data
phased_HO_Euro_1kg_chr19.vcf.gz: genetic data for chromosome 19
probe.pdf: figure showing results of probing
Scripts
interpolate_cM.R: simple script to do linear interpolation on a genetic map.
make_probe_set.sh: makes a "probe" dataset where a chosen region is filled with real data and the rest is filled with fake "IBS inert" genotypes.
write_inert_ped.R: make a fake "IBS-inert" dataset (called by make_probe_set)
analyze_IBDprobe*.R: check IBS probing success

