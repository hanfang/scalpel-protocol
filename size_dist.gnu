reset
n=200	#number of intervals
max=100 #max value
min=-100 #min value
binwidth=(max-min)/n	#interval width

#function used to map a value to the intervals
#hist(x,width)=width*floor(x/width)+width/2.0

set terminal pdf transparent enhanced font 'Helvetica,6'
set output outfile
#set output 'indel_size_dist.pdf'
set xrange [min:max]
set yrange [*:*]

#to put an empty boundary around the
#data inside an autoscaled graph.
set offset graph 0.05,0.05,0.05,0.0
set xtics min,(max-min)/5,max
#set boxwidth width*1.0
#set boxwidth 0.9 absolute

#set style data histogram
#set style data linespoints
#set style histogram rowstacked
#set style histogram clustered gap 1
set style fill solid border -1

#binwidth=5
set boxwidth binwidth*0.9

#function used to map a value to the intervals
bin(x,width)=width*floor(x/binwidth) + binwidth/2.0

#set style line 1 lt 1 lw 3 pt 3
set style fill solid 0.5	#fillstyle
set tics out nomirror
#set title "Exome sample K24510-88962 (known indels -dbSNP 135)"
set xlabel "Size (base pairs)"
set ylabel "Frequency"
set xtic 20

set multiplot layout 1, 1 title ""
#set multiplot layout 1, 3 title "Indel size distribution"
set tmargin 4

#set key box
set key at 100, 8000
show key

# set title "Inherited and de novo indels"
# set label 1 "Scalpel" at -100,400 center

unset logscale y
set table "all.table"
plot infile u (bin($1,binwidth)):(1.0) smooth freq
#plot "all.indel.size" u (bin($1,binwidth)):(1.0) smooth freq
# set table "denovo.table"
# plot "denovo.annovar.filter.ms.noms.frame.exonic.size" u (bin($1,binwidth)):(1.0) smooth freq
unset table
set logscale y
#plot "scalpel.table" smooth freq w boxes title ""
#plot "scalpel.table" smooth freq title "All" lt 9 lw 1 w filledcurves x1
plot "all.table" u 1:($2+1) smooth freq title "Inherited and denovo indels" lt 9 lw 1 w boxes # , \
# "denovo.table" u 1:($2+1) smooth freq title "de novo indels" lt 8 lw 1 w boxes

#plot "scalpel.K8101-49685s.indels.intarget.normal.annovar.size" u (bin($1,binwidth)):(1.0) smooth freq w boxes title ""

#plot "haplotyecaller.K8101-49685s.indels.intarget.normal.annovar.size" u (bin($1,binwidth)):(1.0) smooth freq w boxes title ""

unset multiplot
