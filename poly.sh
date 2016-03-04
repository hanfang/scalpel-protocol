#!/bin/sh


awk -F "\t" '{if($0 ~ /^#/){print $0} else {if( ($7~/LowAltCntAff/ && $7~/HighChi2score/) || $7~/LowCovUnaff/) ) print} }' denovo.twopass.vcf.ms > denovo.twopass.vcf.ms.lq
for i in A C G T; do awk -v j=$i '$0!~/^#/  {  if($15==j) { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }' denovo.twopass.vcf.ms.lq >  poly${i}.VAF.txt ; done
gnuplot44 -e "outfile='homo.vaf.pdf'; infileA='polyA.VAF.txt'; infileC='polyC.VAF.txt'; infileG='polyG.VAF.txt'; infileT='polyT.VAF.txt' " hp.vafdist.gnu

# cat denovo.twopass.vcf.ms inherited.onepass.vcf.ms | grep -v '#' |grep 'yes' > combine.ms.txt
# for i in A C G T; do awk -v j=$i '$0!~/^#/  {  if($15==j) { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }' combine.ms.txt >  poly${i}.VAF.txt ; done 
# gnuplot44 -e "outfile='homo.vaf.pdf'; infileA='polyA.VAF.txt'; infileC='polyC.VAF.txt'; infileG='polyG.VAF.txt'; infileT='polyT.VAF.txt' " hp.vafdist.gnu





# awk '$0!~/^#/ {if($15=="A") { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }'  $input > polyA.VAF.txt
# awk '$0!~/^#/ {if($15=="C") { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }'  $input > polyC.VAF.txt
# awk '$0!~/^#/ {if($15=="G") { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }'  $input > polyG.VAF.txt
# awk '$0!~/^#/ {if($15=="T") { split($12,a,":"); if(a[1]=="0/1" || a[1]=="1/1") split(a[2],b,","); print b[1] "\t" b[2]} }'  $input > polyT.VAF.txt
