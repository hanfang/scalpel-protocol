###----------------------------------------------------------
# README 
### Scalpel resource bundle for variant calling and analysis
###----------------------------------------------------------
### version: v0.5.2
### Authors: Han Fang, Giuseppe Narzisi, Michael C. Schatz
### Date: Sept 4, 2015
###----------------------------------------------------------

#### Introduction:

This resource bundle includes resources for indel variant analysis of whole genome and exome capture sequencing experiments.
Currently we host a few relevant files for using scalpel (http://scalpel.sourceforge.net).

#### Resources:

- msdetector: a micro-satellite detector for a given vcf file (marks variants within STR regions)
- SeqCap_EZ_Exome_v3_primary.scalpel.bed: exonic targetted bed file derived from the bed file of SeqCap EZ Human Exome Library v3.0 (http://www.nimblegen.com/downloads/annotation/ez_exome_v3/SeqCapEZ_Exome_v3.0_Design_Annotation_files.zip)
- clinvar_main.bed: a bed file of the variants in the ClinVar main database June 2015 
- *gnu: several gnuplot files for plotting
- filter_vcf.sh: a shell wrapped for basic filtering of the scalpel vcf files
- run_protocol.sh: a shell script documenting the main step in the protocol manuscript


#### Cite:
- Fang H, Grabowska EA, Arora K, Vacic V, Zody M, Iossifov I, O’Rawe JA, Y Wu, Jimenez-Barron LT, Rosenbaum J, Ronemus M, Lee Y, Wang Z, Gholson J. Lyon, Michael Wigler, Michael C. Schatz, Giuseppe Narzisi, "De novo and somatic indel variant analysis of whole genome and exome capture sequencing experiments with Scalpel", (In preparation)
- Fang H, Wu Y, Narzisi G, O’Rawe JA, Jimenez Barron LT, Rosenbaum J, Ronemus M, Iossifov I, Schatz MC*, Lyon GJ* (2014). " Reducing INDEL calling errors in whole-genome and exome sequencing data". Genome Medicine (2014) doi:10.1186/s13073-014-0089-z 
- Narzisi G, O’Rawe JA, Iossifov I, Fang H, Lee Y, Wang Z, Wu Y, Lyon GJ, Wigler M, Schatz MC. "Accurate detection of de novo and transmitted INDELs within exome-capture data using micro-assembly". Nature Methods (2014) doi:10.1038/nmeth.3069. 