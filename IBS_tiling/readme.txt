Scripts to generate IBS tiling results for "Attacks on genetic privacy via uploads to genealogical databases" by Michael D. Edge and Graham Coop

These scripts require the ability to run R, java, germline, EIGENSTRAT, and plink from the command line.

Software included:

beagle.16May19.351.jar - Executable for Beagle phasing software.
refined-ibd.16May19.jar - Executable for RefinedIBD program used for calling IBD / IBS segments from data phased by Beagle.

Scripts:

merge_HO_and1kg.sh - merges genetic data from 1000 Genomes Project and Reich lab
par.PED.EIGENSTRAT - used to call EIGENSTRAT as part of processing genetic data from Reich lab
phase_ibs_HO_and_1kg.sh - calls Beagle to phase data and calls RefinedIBD to identify IBD / IBS segments
call_ibs_germline.sh - calls germline to identify IBD / IBS segments using germline
IBStiling*.R - R script that identifies segments "tiled" according to IBD/IBS results and make figures. These are mostly the same, with the important difference being the IBD segments used as input.

The recombination map we used is provided in plink.GRCh37.map; it was originally obtained from Brian Browning's Beagle website.
