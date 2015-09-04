#!/bin/sh
set -e
###----------------------------------------------------------
### Scalpel resource bundle for variant calling and analysis
###----------------------------------------------------------
### version: v0.5.2
### Authors: Han Fang, Giuseppe Narzisi, Michael C. Schatz
### Date: Sept 4, 2015
###----------------------------------------------------------

##--------------------------------
## download and set up files/tools
##--------------------------------
## 1|Download the example sequencing reads of the Hapmap quad family from the Illumina Platinum Genome project (*_1*fastq.gz and *_2*fastq.gz denote paired end reads):
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_2.fastq.gz

## 2|Download the human reference genome hg19:
wget --no-check http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
wget --no-check http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

## 3|Convert the *.2bit genome to *.fa format and index it with bwa (Note you can also download the fasta file directly, although may take much longer):
chmod +x twoBitToFa ; ./twoBitToFa hg19.2bit hg19.fa
bwa index hg19.fa

##--------------------------------
## Align the NGS reads to the genome
##--------------------------------
## 4|Align reads to reference for each sample separately with bwa mem:
bwa mem -t 10 -R '@RG\tID:NA12877\tSM:NA12877' hg19.fa ERR194146_1.fastq.gz ERR194146_2.fastq.gz | samtools view -h -S -b > NA12877.bam
bwa mem -t 10 -R '@RG\tID:NA12878\tSM:NA12878' hg19.fa ERR194147_1.fastq.gz ERR194147_2.fastq.gz | samtools view -h -S -b > NA12878.bam
bwa mem -t 10 -R '@RG\tID:NA12881\tSM:NA12881' hg19.fa ERR324432_1.fastq.gz ERR324432_2.fastq.gz | samtools view -h -S -b > NA12881.bam
bwa mem -t 10 -R '@RG\tID:NA12882\tSM:NA12882' hg19.fa ERR194151_1.fastq.gz ERR194151_2.fastq.gz | samtools view -h -S -b > NA12882.bam

## 5|Sort the bam files by chromosome coordinates with samtools:
samtools sort -m 4G NA12877.bam NA12877.sort
samtools sort -m 4G NA12878.bam NA12878.sort
samtools sort -m 4G NA12881.bam NA12881.sort
samtools sort -m 4G NA12882.bam NA12882.sort
rm -f NA12877.bam NA12878.bam NA12881.bam NA12882.bam 

## 6|Mark duplicated reads within the alignment with picard tools:
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12877.sort.bam OUTPUT=NA12877.sort.markdup.bam METRICS_FILE=NA12877.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12878.sort.bam OUTPUT=NA12878.sort.markdup.bam METRICS_FILE=NA12878.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12881.sort.bam OUTPUT=NA12881.sort.markdup.bam METRICS_FILE=NA12881.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12882.sort.bam OUTPUT=NA12882.sort.markdup.bam METRICS_FILE=NA12882.sort.metric
rm -f NA12877.sort.bam NA12878.sort.bam NA12881.sort.bam NA12882.sort.bam

## 7|Perform a basic quality control of the alignment files with samtools:
samtools flagstat NA12877.sort.markdup.bam > NA12877.sort.markdup.bam.simplestats
samtools flagstat NA12878.sort.markdup.bam > NA12878.sort.markdup.bam.simplestats
samtools flagstat NA12881.sort.markdup.bam > NA12881.sort.markdup.bam.simplestats
samtools flagstat NA12882.sort.markdup.bam > NA12882.sort.markdup.bam.simplestats

##--------------------------------
## Perform indel variant calling and downstream filtering
##--------------------------------
## 8|Run Scalpel in the “de novo” mode to perform multi-sample calling for a family. In this example, we use NA12882 as the affected individual. The NA12881 is the unaffected individual accordingly:
scalpel-discovery --denovo --dad NA12877.sort.markdup.bam --mom NA12878.sort.markdup.bam --aff NA12882.sort.markdup.bam --sib NA12881.sort.markdup.bam --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --numprocs 10 --two-pass

## 9|Export the inherited and denovo mutations from the Scalpel database (in target only):
scalpel-export --denovo --db outdir/main/inherited.db  --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget --min-alt-count-affected 10 --max-chi2-score 10.8  > inherited.onepass.vcf
scalpel-export --denovo --db outdir/twopass/denovos.db --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget --min-alt-count-affected 10 --max-chi2-score 10.8 --min-coverage-unaffected 20 > denovo.twopass.vcf

## 10|Identify and mark indels within STR regions using ms-detector:
sh ./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i inherited.onepass.vcf > inherited.onepass.vcf.ms
sh ./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i denovo.twopass.vcf    > denovo.twopass.vcf.ms

## 11|Save indels within and outside STR regions into different vcf files:
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="yes") print} }' inherited.onepass.vcf.ms | cut -f1-13 > inherited.onepass.vcf.ms.in
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="no") print} }' inherited.onepass.vcf.ms  | cut -f1-13  > inherited.onepass.vcf.ms.out
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="yes") print} }' denovo.twopass.vcf.ms    | cut -f1-13 > denovo.twopass.vcf.ms.in
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="no") print} }' denovo.twopass.vcf.ms     | cut -f1-13 > denovo.twopass.vcf.ms.out

## 12|Filter out false positive calls by adjusting coverage and/or chi-squared thresholds for your data:
awk -F "\t" '{if($0 ~ /^#/){print $0} else {if(! ($7~/LowAltCntAff/ && $7~/HighChi2score/) ) print} }' inherited.onepass.vcf.ms.out > inherited.onepass.vcf.ms.out.hq
awk -F "\t" '{if($0 ~ /^#/){print $0} else {if(! ($7~/LowAltCntAff/ || $7~/HighChi2score/ || $7~/LowCovUnaff/) ) print} }' denovo.twopass.vcf.ms.out > denovo.twopass.vcf.ms.out.hq

## 13|(Optional) Further customize your variant call set with a python script 
python denovo-multi-filter.py -i denovo.twopass.vcf.ms.out -f NA12877 -m NA12878 -a NA12882 -u NA12881 -aac 10 -chi 10.8  -pc 20 -o denovo.twopass.vcf.ms.out.filter

## 14|(Optional) Extract a subset of indels based on other annotation using bedtools:
bedtools intersect -wa -u -a inherited.onepass.vcf.ms.out.hq -b clinvar_main.bed > inherited.onepass.vcf.ms.out.hq.clinvar

## 15|Summarize indel calls with histogram of mutations by size.
grep -v "#" inherited.onepass.vcf.ms.out.hq denovo.twopass.vcf.ms.out.hq | awk '{print length($5)-length($4)}' > all.indel.size.txt
gnuplot44 -e "outfile='indel_size_dist.pdf'; infile='all.indel.size.txt'" size_dist.gnu

## 16|Summarize homopolymer indels calls with histogram of mutations by VAF.
cat denovo.twopass.vcf.ms inherited.onepass.vcf.ms | grep -v '#' |grep 'yes' | awk -F "\t" '{if( ($7~/LowAltCntAff/ && $7~/HighChi2score/) || $7~/LowCovUnaff/ ) print}' > combine.ms.txt
for i in A C G T; do awk -v j=$i '$0!~/^#/  {  if($15==j) { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }' combine.ms.txt >  poly${i}.VAF.txt ; done
gnuplot44 -e "outfile='homo.vaf.pdf'; infileA='polyA.VAF.txt'; infileC='polyC.VAF.txt'; infileG='polyG.VAF.txt'; infileT='polyT.VAF.txt' " hp.vafdist.gnu

## 17|Summarize inherited indels with variant allele frequencies (VAF %): 
awk -F'\t' '$0!~/^#/ {split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]}' inherited.onepass.vcf.ms.out > inherited.onepass.vcf.ms.out.vaf
awk -F'\t' '$0!~/^#/ {split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]}' inherited.onepass.vcf.ms.out.hq > inherited.onepass.vcf.ms.out.hq.vaf
gnuplot44 -e "outfile='inherited.VAFdist.pdf'; infileAll='inherited.onepass.vcf.ms.out.vaf'; infileHq='inherited.onepass.vcf.ms.out.hq.vaf'" vafdistplot.inherited.qual.gnu

## 18|Determine the number of indels remained after each step of the filtering:
for i in *.vcf.* ; do echo $i; grep -v "#" $i | wc -l ;done

## 19|split the multisample vcf to individual files
for file in *.hq; do bgzip -c $file > $file.gz ; tabix -p vcf $file.gz; done
for file in *.hq.gz; do for sample in `bcftools view -h $file | grep "^#CHROM" | cut -f10-`; do bcftools view -c1 -Oz -s $sample -o ${file/.gz*/.$sample.vcf.gz} ${file}; gunzip ${file/.gz*/.$sample.vcf.gz}; done;done

##--------------------------------
## Annotation and visualization of the indel calls 
##--------------------------------
## 20|Prepare and create the input format required by ANNOVAR: 
# annovar=/path-to-annovar/
$annovar/convert2annovar.pl -format vcf4 inherited.onepass.vcf.ms.out.hq.NA12882.vcf > inherited.onepass.vcf.ms.out.hq.NA12882.vcf.avinput
$annovar/convert2annovar.pl -format vcf4 denovo.twopass.vcf.ms.out.hq.NA12882.vcf    > denovo.twopass.vcf.ms.out.hq.NA12882.vcf.avinput

## 21|Annotate and intersect indels with gene regions using ANNOVAR
$annovar/annotate_variation.pl -buildver hg19 inherited.onepass.vcf.ms.out.hq.NA12882.vcf.avinput $annovar/humandb
$annovar/annotate_variation.pl -buildver hg19 denovo.twopass.vcf.ms.out.hq.NA12882.vcf.avinput $annovar/humandb

## 22|Summarize coding region indels by size in R.
cat inherited.onepass.vcf.ms.out.hq.NA12882.vcf.avinput.exonic_variant_function | egrep -v 'unknown|stopgain' | cut -f 2,7,8 | cut -d " " -f 2 | awk '{if($2=="-") print $1"\t"length($3);else if ($3=="-") print $1"\t"length($2)}' > type_and_size.txt
Rscript coding_indel_size.R
# % R
# > indel=read.table("type_and_size.txt", header=FALSE)
# > colnames(indel)= c("type","size")
# > indel_30=indel[indel[,2]<=30,]
# > indel.table <- table(indel_30$type,factor(indel_30$size,lev=1:30)  )
# > pdf('indelsize_by_type.pdf', width=12, height=8)
# > barplot(indel.table, main="indel distribution within coding sequence (CDS)", xlab="", col=c("green","red"), legend = rownames(indel.table))
# > dev.off()

## 23|Retrieve frame-shift mutations, which are potentially loss-of-function
awk '{if($2=="frameshift") print}' inherited.onepass.vcf.ms.out.hq.NA12882.vcf.avinput.exonic_variant_function > inherited.onepass.vcf.ms.out.hq.NA12882.vcf.avinput.exonic_variant_function.frameshift
awk '{if($2=="frameshift") print}' denovo.twopass.vcf.ms.out.hq.NA12882.vcf.avinput.exonic_variant_function > denovo.twopass.vcf.ms.out.hq.NA12882.vcf.avinput.exonic_variant_function.frameshift