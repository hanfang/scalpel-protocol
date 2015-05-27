#!/bin/sh

################################
## Software setup
################################

## Download some relevant files that are host on github, including msdetector.
# git clone https://github.com/hanfang/scalpel-protocol.git; cd msdetector; make; chmod +x *; cd ../

## Download tools that are required for the indel analysis
## bwa setup
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2
tar jxf bwa-0.7.12.tar.bz2; cd bwa-0.7.12; make; cd ../ ; export PATH=./bwa-0.7.12:$PATH

## samtools setup
wget http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2
tar jxf samtools-1.2.tar.bz2; cd samtools-1.2; make; cd ../ ; export PATH=./samtools-1.2:$PATH

## picard setup
wget https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip --no-check-certificate ; unzip picard-tools-1.130.zip ; export PATH=./picard-tools-1.130:$PATH

## scalpel setup
wget http://sourceforge.net/projects/scalpel/files/scalpel-0.4.1.tar.gz --no-check-certificate ; tar zxvf scalpel-0.4.1.tar.gz
cd scalpel-0.4.1; make; cd ../; export PATH=./scalpel-0.4.1:$PATH

## bedtools setup
wget https://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz --no-check-certificate ; tar zxvf bedtools-2.23.0.tar.gz
cd bedtools2; make; cd ../ ; export PATH=./bedtools-2.23.0/bin:$PATH

## annovar setup
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz --no-check-certificate ; tar zxvf annovar.latest.tar.gz ; # export PATH=./annovar:$PATH

################################
## PROCEDURE
################################

## 1.Downloading sequencing data; Example data (Illumina Platinum Genome)

## Download example fastq files - Hapmap quad family; *_1*fastq.gz and *_2*fastq.gz denote paired end reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194151/ERR194151_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR324/ERR324432/ERR324432_2.fastq.gz

## Reference data set up - download hg19 reference genome (2bit) from UCSC ftp and the twoBitToFa coverter
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

## Convert the 2bit file to fa format
chmod +x twoBitToFa ; ./twoBitToFa hg19.2bit hg19.fa


## 2.Align reads to reference for each sample separately (align, sort, index, etc.)

## Index reference genome
bwa index hg19.fa

## align the reads with bwa mem, convert the sam file to bam file with samtools
bwa mem -t 10 hg19.fa ERR194146_1.fastq.gz ERR194146_2.fastq.gz | samtools view -h -S -b > NA12877.bam
bwa mem -t 10 hg19.fa ERR194147_1.fastq.gz ERR194147_2.fastq.gz | samtools view -h -S -b > NA12878.bam
bwa mem -t 10 hg19.fa ERR324432_1.fastq.gz ERR324432_2.fastq.gz | samtools view -h -S -b > NA12881.bam
bwa mem -t 10 hg19.fa ERR194151_1.fastq.gz ERR194151_2.fastq.gz | samtools view -h -S -b > NA12882.bam

## sort the bam files
samtools sort -m 4G NA12877.bam NA12877.sort
samtools sort -m 4G NA12878.bam NA12878.sort
samtools sort -m 4G NA12881.bam NA12881.sort
samtools sort -m 4G NA12882.bam NA12882.sort
rm -f NA12877.bam NA12878.bam NA12881.bam NA12882.bam 

## mark duplicates with picard tools
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12877.sort.bam OUTPUT=NA12877.sort.markdup.bam METRICS_FILE=NA12877.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12878.sort.bam OUTPUT=NA12878.sort.markdup.bam METRICS_FILE=NA12878.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12881.sort.bam OUTPUT=NA12881.sort.markdup.bam METRICS_FILE=NA12881.sort.metric
java -jar -Xmx10g picard MarkDuplicates INPUT=NA12882.sort.bam OUTPUT=NA12882.sort.markdup.bam METRICS_FILE=NA12882.sort.metric
rm -f NA12877.sort.bam NA12878.sort.bam NA12881.sort.bam NA12882.sort.bam

## 3.Run Scalpel in “de novo” mode to perform multi-sample calling for a family

## Use scalpel to call indels
scalpel --denovo --dad NA12877.sort.markdup.bam --mom NA12878.sort.markdup.bam --aff NA12882.sort.markdup.bam –sib NA12881.sort.markdup.bam --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --numprocs 10 --two-pass

## 4.Export mutations (--export)

## Export inherited indels - in target only
scalpel --export --db outdir/main/inherited.db --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget > inherited.onepass.vcf

## Export denovo indels - in target only
scalpel --export --db outdir/twopass/denovos.db --bed SeqCap_EZ_Exome_v3_primary.scalpel.bed --ref hg19.fa --intarget > denovo.twopass.vcf

## 5.Summarize data and generate QC Plots (indel size, ref-alt cov, etc)

## 6.Critical issue: much higher sequencing biases in GC-extreme regions. Detect indels within STRs (MSdetector) and near identical repeats. 

## identify indels within STR regions with ms-detector
./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i inherited.onepass.vcf > inherited.onepass.vcf.ms
./msdetector/msdetector.sh -r 50 -d 2 -g hg19.fa -i denovo.twopass.vcf    > denovo.twopass.vcf.ms

## Save indels within and outside STR regions into different vcf files
cat inherited.onepass.vcf.ms | awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($13=="yes") print} }' | cut -d\t -f1-10 > inherited.onepass.vcf.ms.in
cat inherited.onepass.vcf.ms | awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($13=="no") print} }' | cut -d\t -f1-10  > inherited.onepass.vcf.ms.out
cat denovo.twopass.vcf.ms    | awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($13=="yes") print} }' | cut -d\t -f1-10 > denovo.twopass.vcf.ms.in
cat denovo.twopass.vcf.ms    | awk -F "\t" '{if($0 ~ /^#/){print $0} else{if($13=="no") print} }' | cut -d\t -f1-10 > denovo.twopass.vcf.ms.out


## 7.Critical Issue: Remove false-positives by adjusting coverage and/or chi-squared thresholds for your data.

## filter out low-quality indels
cat inherited.onepass.vcf.ms.out | awk -F "[\t;=,/]" '{ if($0 ~ /^#/) {print $0}  else  { if   ($11 >=10 || $19 <=10.8) print } } ' > inherited.onepass.vcf.ms.out.hq
cat denovo.twopass.vcf.ms.out    | awk -F "[\t;=,/]" '{ if($0 ~ /^#/) {print $0}  else  { if ( ($11 >=10 || $19 <=10.8) && ($26 >=20 && $27>=20) ) print } } ' > denovo.twopass.vcf.ms.out.hq


## 8.Interpreting results: intersect indels with gene regions using annovar, retrieve population frequencies of indels from 1000g and exome sequencing project database

## annotate the indels with annovar
./annovar/convert2annovar.pl -format vcf4 inherited.onepass.vcf.ms.out.hq > inherited.onepass.vcf.ms.out.hq.avinput
./annovar/annotate_variation.pl -buildver hg19 inherited.onepass.vcf.ms.out.hq.avinput ./annovar/humandb

./annovar/convert2annovar.pl -format vcf4 denovo.twopass.vcf.ms.out.hq > denovo.twopass.vcf.ms.out.hq.avinput
./annovar/annotate_variation.pl -buildver hg19 denovo.twopass.vcf.ms.out.hq.avinput ./annovar/humandb
