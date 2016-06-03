#!/bin/sh
set -e

##----------------------------------------------------------
## Scalpel resource bundle for variant calling and analysis
##----------------------------------------------------------
## version: v0.5.3
## Authors: Han Fang, Giuseppe Narzisi, Michael C. Schatz
## Date: March 4, 2016
##----------------------------------------------------------

##--------------------------------
## define your variables here
##--------------------------------
## These variables will serve as prefixes for files in the analysis. If you wish to analyze your own data, please change the variables here.
father=NA12877
mother=NA12878
proband=NA12882
sibling=NA12881

##--------------------------------
## download and set up files/tools
##--------------------------------
## 1| Download the example sequencing reads of the Hapmap quad family from the Illumina Platinum Genome project (*_1*fastq.gz and *_2*fastq.gz denote paired end reads):
## You can skip this if you decide to use your own fastq files, rather than the example data.
echo "[status] 1| Download the example sequencing reads"
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_2.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_1.fastq.gz
wget --no-check ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_2.fastq.gz

## 2| Download the human reference genome hg19:
## You can skip this if you already have a reference genome downloaded.
echo "[status] 2| Download the human reference genome hg19"
wget --no-check http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
wget --no-check http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

## 3| Convert the *.2bit genome to *.fa format and index it with bwa (Note you can also download the fasta file directly, although may take much longer).
## You can skip this if you already have a reference genome indexed.
echo "[status] 3| Convert the *.2bit genome to *.fa format and index it"
chmod +x twoBitToFa ; ./twoBitToFa hg19.2bit hg19.fa
bwa index hg19.fa

##--------------------------------
## Align the NGS reads to the genome
##--------------------------------
## 4| Align reads to reference for each sample separately with bwa mem:
echo "[status] 4| Align reads to reference"
bwa mem -t 10 -R '@RG\tID:'${father}'\tSM:'${father} hg19.fa ERR194146_1.fastq.gz ERR194146_2.fastq.gz | samtools view -h -S -b > ${father}.bam
bwa mem -t 10 -R '@RG\tID:'${mother}'\tSM:'${mother} hg19.fa ERR194147_1.fastq.gz ERR194147_2.fastq.gz | samtools view -h -S -b > ${mother}.bam
bwa mem -t 10 -R '@RG\tID:'${sibling}'\tSM:'${sibling} hg19.fa ERR324432_1.fastq.gz ERR324432_2.fastq.gz | samtools view -h -S -b > ${sibling}.bam
bwa mem -t 10 -R '@RG\tID:'${proband}'\tSM:'${proband} hg19.fa ERR194151_1.fastq.gz ERR194151_2.fastq.gz | samtools view -h -S -b > ${proband}.bam

## 5| Sort the bam files by chromosome coordinates with samtools:
echo "[status] 5| Sort the bam files"
samtools sort -m 4G -o ${father}.sort.bam ${father}.bam
samtools sort -m 4G -o ${mother}.sort.bam ${mother}.bam
samtools sort -m 4G -o ${sibling}.sort.bam ${sibling}.bam
samtools sort -m 4G -o ${proband}.sort.bam ${proband}.bam
rm -f ${father}.bam ${mother}.bam ${sibling}.bam ${proband}.bam

## 6| Mark duplicated reads within the alignment with picard tools:
echo "[status] 6| Mark duplicated reads"
java -jar -Xmx10g picard.jar MarkDuplicates INPUT=${father}.sort.bam OUTPUT=${father}.sort.markdup.bam METRICS_FILE=${father}.sort.metric
java -jar -Xmx10g picard.jar MarkDuplicates INPUT=${mother}.sort.bam OUTPUT=${mother}.sort.markdup.bam METRICS_FILE=${mother}.sort.metric
java -jar -Xmx10g picard.jar MarkDuplicates INPUT=${sibling}.sort.bam OUTPUT=${sibling}.sort.markdup.bam METRICS_FILE=${sibling}.sort.metric
java -jar -Xmx10g picard.jar MarkDuplicates INPUT=${proband}.sort.bam OUTPUT=${proband}.sort.markdup.bam METRICS_FILE=${proband}.sort.metric
rm -f ${father}.sort.bam ${mother}.sort.bam ${sibling}.sort.bam ${proband}.sort.bam

## 7| Perform a basic quality control of the alignment files with samtools:
echo "[status] 7| Perform a basic quality control"
samtools flagstat ${father}.sort.markdup.bam > ${father}.sort.markdup.bam.simplestats
samtools flagstat ${mother}.sort.markdup.bam > ${mother}.sort.markdup.bam.simplestats
samtools flagstat ${sibling}.sort.markdup.bam > ${sibling}.sort.markdup.bam.simplestats
samtools flagstat ${proband}.sort.markdup.bam > ${proband}.sort.markdup.bam.simplestats

##--------------------------------
## Perform indel variant calling and downstream filtering
##--------------------------------
## 8| Run Scalpel in the “de novo” mode to perform multi-sample calling for a family. In this example, we use NA12882 as the proband (affected individual). The NA12881 is the unaffected sibling accordingly:
echo "[status] 8| Run Scalpel in the “de novo” mode to perform multi-sample calling for a family"
scalpel-discovery --denovo --dad ${father}.sort.markdup.bam --mom ${mother}.sort.markdup.bam --aff ${proband}.sort.markdup.bam --sib ${sibling}.sort.markdup.bam --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --numprocs 10 --two-pass

## 9| Export the inherited and denovo mutations from the Scalpel database (in target only):
echo "[status] 9| Export the inherited and denovo mutations from the Scalpel database"
scalpel-export --denovo --db outdir/main/inherited.db  --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget --min-alt-count-affected 10 --max-chi2-score 10.8  > inherited.onepass.vcf
scalpel-export --denovo --db outdir/twopass/denovos.db --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget --min-alt-count-affected 10 --max-chi2-score 10.8 --min-coverage-unaffected 20 > denovo.twopass.vcf

## 10| Identify and mark indels within STR regions using ms-detector:
echo "[status] 10| Identify and mark indels within STR regions"
sh ./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i inherited.onepass.vcf > inherited.onepass.vcf.ms
sh ./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i denovo.twopass.vcf    > denovo.twopass.vcf.ms

## 11| Save indels within and outside STR regions into different vcf files:
echo "[status] 11| Save indels within and outside STR regions into different vcf files"
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="yes") print} }' inherited.onepass.vcf.ms | cut -f1-13 > inherited.onepass.vcf.ms.in
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="no") print} }' inherited.onepass.vcf.ms  | cut -f1-13  > inherited.onepass.vcf.ms.out
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="yes") print} }' denovo.twopass.vcf.ms    | cut -f1-13 > denovo.twopass.vcf.ms.in
awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($16=="no") print} }' denovo.twopass.vcf.ms     | cut -f1-13 > denovo.twopass.vcf.ms.out

## 12| Filter out false positive calls by adjusting coverage and/or chi-squared thresholds for your data:
echo "[status] 12| Filter out false positive calls by adjusting coverage and/or chi-squared thresholds"
awk -F "\t" '{if($0 ~ /^#/){print $0} else {if(! ($7~/LowAltCntAff/ && $7~/HighChi2score/) ) print} }' inherited.onepass.vcf.ms.out > inherited.onepass.vcf.ms.out.hq
awk -F "\t" '{if($0 ~ /^#/){print $0} else {if(! ($7~/LowAltCntAff/ || $7~/HighChi2score/ || $7~/LowCovUnaff/) ) print} }' denovo.twopass.vcf.ms.out > denovo.twopass.vcf.ms.out.hq

## 13| (Optional) Further customize your variant call set with a python script
echo "[status] 13| (Optional) Further customize your variant call set"
python denovo-multi-filter.py -i denovo.twopass.vcf.ms.out -f ${father} -m ${mother} -a ${proband} -u ${sibling} -aac 10 -chi 10.8  -pc 20 -o denovo.twopass.vcf.ms.out.filter

## 14| (Optional) Extract a subset of indels based on other annotation using bedtools:
echo "[status] 14| (Optional) Extract a subset of indels based on other annotation"
bedtools intersect -wa -u -a inherited.onepass.vcf.ms.out.hq -b clinvar_main.bed > inherited.onepass.vcf.ms.out.hq.clinvar

## 15| Summarize indel calls with histogram of mutations by size.
echo "[status] 15| Summarize indel calls with histogram of mutations by size"
grep -v "#" inherited.onepass.vcf.ms.out.hq denovo.twopass.vcf.ms.out.hq | awk '{print length($5)-length($4)}' > all.indel.size.txt
gnuplot44 -e "outfile='indel_size_dist.pdf'; infile='all.indel.size.txt'" size_dist.gnu

## 16| Summarize homopolymer indels calls with histogram of mutations by VAF.
echo "[status] 16| Summarize homopolymer indels calls with histogram of mutations by VAF"
cat denovo.twopass.vcf.ms inherited.onepass.vcf.ms | grep -v '#' |grep 'yes' | awk -F "\t" '{if( ($7~/LowAltCntAff/ && $7~/HighChi2score/) || $7~/LowCovUnaff/ ) print}' > combine.ms.txt
for i in A C G T; do awk -v j=$i '$0!~/^#/  {  if($15==j) { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }' combine.ms.txt >  poly${i}.VAF.txt ; done
gnuplot44 -e "outfile='homo.vaf.pdf'; infileA='polyA.VAF.txt'; infileC='polyC.VAF.txt'; infileG='polyG.VAF.txt'; infileT='polyT.VAF.txt' " hp.vafdist.gnu

## 17| Summarize inherited indels with variant allele frequencies (VAF %): 
echo "[status] 17| Summarize inherited indels with variant allele frequencies (VAF %)"
awk -F'\t' '$0!~/^#/ {split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]}' inherited.onepass.vcf.ms.out > inherited.onepass.vcf.ms.out.vaf
awk -F'\t' '$0!~/^#/ {split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]}' inherited.onepass.vcf.ms.out.hq > inherited.onepass.vcf.ms.out.hq.vaf
gnuplot44 -e "outfile='inherited.VAFdist.pdf'; infileAll='inherited.onepass.vcf.ms.out.vaf'; infileHq='inherited.onepass.vcf.ms.out.hq.vaf'" vafdistplot.inherited.qual.gnu

## 18| Determine the number of indels remained after each step of the filtering:
echo "[status] 18| Determine the number of indels remained after each step of the filtering"
for i in *.vcf.* ; do echo $i; grep -v "#" $i | wc -l ;done > indel.count.txt

## 19| Split the multisample vcf to an individual file for the proband
echo "[status] 19| Split the multisample vcf to an individual file for the proband"
for file in *.hq; do bgzip -c $file > $file.gz ; tabix -p vcf $file.gz; done
for file in *.hq.gz; do bcftools view -c1 -Ov -s ${proband} -o ${file/.gz*/.${proband}.vcf} ${file}; done

## 20| Filter the single vcf files based on Chi-Square score and allele coverage
echo "[status] 20| Filter the single vcf files based on Chi-Square score and allele coverage"
python single-vcf-filter.py -i inherited.onepass.vcf.ms.out.hq.${proband}.vcf -mc 10 -chi 10.8 -o inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf
python single-vcf-filter.py -i denovo.twopass.vcf.ms.out.hq.${proband}.vcf -mc 10 -chi 10.8 -o denovo.twopass.vcf.ms.out.hq.${proband}.filter.vcf

##--------------------------------
## Annotation and visualization of the indel calls 
##--------------------------------
## 21| Prepare and create the input format required by ANNOVAR: 
echo "[status] 21| Prepare and create the input format required by ANNOVAR"
annovar/convert2annovar.pl -format vcf4 inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf > inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf.avinput
annovar/convert2annovar.pl -format vcf4 denovo.twopass.vcf.ms.out.hq.${proband}.filter.vcf    > denovo.twopass.vcf.ms.out.hq.${proband}.filter.vcf.avinput

## 22| Annotate and intersect indels with gene regions using ANNOVAR
echo "[status] 22| Annotate and intersect indels with gene regions using ANNOVAR"
annovar/annotate_variation.pl -buildver hg19 inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf.avinput $annovar/humandb
annovar/annotate_variation.pl -buildver hg19 denovo.twopass.vcf.ms.out.hq.${proband}.filter.vcf.avinput $annovar/humandb

## 23| Summarize coding region indels by size in R.
echo "[status] 23| Summarize coding region indels by size in R"
cat inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf.avinput.exonic_variant_function | egrep -v 'unknown|stopgain' | cut -f 2,7,8 | cut -d " " -f 2 | awk '{if($2=="-") print $1"\t"length($3);else if ($3=="-") print $1"\t"length($2)}' > type_and_size.txt
Rscript coding_indel_size.R
# % R
# > indel=read.table("type_and_size.txt", header=FALSE)
# > colnames(indel)= c("type","size")
# > indel_30=indel[indel[,2]<=30,]
# > indel.table <- table(indel_30$type,factor(indel_30$size,lev=1:30)  )
# > pdf('indelsize_by_type.pdf', width=16, height=7)
# > mar.default <- c(5,4,4,2) + 0.1
# > par(mar = mar.default + c(0, 4, 0, 0))
# > barplot(indel.table, main="indel distribution within coding sequence (CDS)", xlab="indel size", ylab="number of indels", col=c("green","red"), cex.axis=2, cex.names=2 , cex.lab = 2 , cex.main=2, cex.sub=2 )
# > legend('topright',rownames(indel.table), fil=c('green', 'red'),  bty='n', cex=2)

## 24| Filter the indels based on population allele frequencies:
echo "[status] 24| Filter the indels based on population allele frequencies"
annovar/annotate_variation.pl -filter -out inherited.onepass.vcf.ms.out.hq.${proband}.filter -dbtype popfreq_max_20150413 -build hg19 inherited.onepass.vcf.ms.out.hq.${proband}.filter.vcf.avinput $annovar/humandb/
annovar/annotate_variation.pl -filter -out denovo.twopass.vcf.ms.out.hq.${proband}.filter -dbtype popfreq_max_20150413 -build hg19 denovo.twopass.vcf.ms.out.hq.${proband}.filter.vcf.avinput $annovar/humandb/

## 25| Annotate novel indels that were not reported by population database before (1000G, ESP6500, ExAC, CG46)
echo "[status] 25| Annotate novel indels that were not reported by population database before"
annovar/annotate_variation.pl -buildver hg19 inherited.onepass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered $annovar/humandb
annovar/annotate_variation.pl -buildver hg19 denovo.twopass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered $annovar/humandb

## 26| Retrieve novel frame-shift mutations, which are potentially loss-of-function
echo "[status] 26| Retrieve novel frame-shift mutations, which are potentially loss-of-function"
awk '{if($2=="frameshift") print}' inherited.onepass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered.exonic_variant_function > inherited.onepass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered.exonic_variant_function.fs.txt
awk '{if($2=="frameshift") print}' denovo.twopass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered.exonic_variant_function > denovo.twopass.vcf.ms.out.hq.${proband}.filter.hg19_popfreq_max_20150413_filtered.exonic_variant_function.fs.txt
