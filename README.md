# README 
#### Scalpel resource bundle for variant calling and analysis
#### version: v0.5.3
#### Authors: Han Fang, Giuseppe Narzisi, Michael C. Schatz
#### Date: June 2, 2016
###----------------------------------------------------------

#### Introduction:

This resource bundle includes resources for indel variant analysis of whole genome and exome capture sequencing experiments.
Currently we host some relevant scripts and files for running Scalpel and analyzing its outputs (http://scalpel.sourceforge.net).

#### Resources:
- Docker image: https://hub.docker.com/r/hanfang/scalpel/
  - ` docker pull hanfang/scalpel `
- Scalpel is also deployed on FireCloud: https://software.broadinstitute.org/firecloud/

- This tar ball contains:
  - run_protocol_0.53.sh: a master shell script documenting the main step in the protocol manuscript
  - msdetector: a micro-satellite detector for marking variants within STR regions in a given vcf file
  - SeqCap_EZ_Exome_v3_primary.scalpel.bed: a bed file contains exonic regions derived from the bed file of SeqCap EZ Human Exome Library v3.0 (http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip)
  - clinvar_main.bed: a bed file of the variants in the ClinVar main database June 2015 
  - *.gnu: several gnuplot scripts for plotting
  - denovo-multi-filter.py: a python script for filtering Scalpel multi sample vcf (Python 2-3 compatible)
  - single-vcf-filter.py: a python script for filtering Scalpel single sample vcf (Python 2-3 compatible)
  - coding_indel_size.R: a R scipt for plotting coding region indel distribution


#### If you use Scalpel, please cite:
- Fang H, Grabowska EA, Arora K, Vacic V, Zody M, Iossifov I, O’Rawe JA, Y Wu, Jimenez-Barron LT, Rosenbaum J, Ronemus M, Lee Y, Wang Z, Dikoglu E, Jobanputra V, Lyon GJ, Wigler M, Schatz MC, Narzisi G*, "Indel variant analysis of short-read sequencing data with Scalpel". Preprint in bioRxiv (2016) doi: dx.doi.org/10.1101/028050
- Fang H, Wu Y, Narzisi G, O’Rawe JA, Jimenez Barron LT, Rosenbaum J, Ronemus M, Iossifov I, Schatz MC*, Lyon GJ* (2014). " Reducing INDEL calling errors in whole-genome and exome sequencing data". Genome Medicine (2014) doi:10.1186/s13073-014-0089-z 
- Narzisi G, O’Rawe JA, Iossifov I, Fang H, Lee Y, Wang Z, Wu Y, Lyon GJ, Wigler M, Schatz MC. "Accurate detection of de novo and transmitted INDELs within exome-capture data using micro-assembly". Nature Methods (2014) doi:10.1038/nmeth.3069.
