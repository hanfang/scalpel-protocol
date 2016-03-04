reset
n=100	#number of intervals
max=100 #max value
min=0 #min value
binwidth=(max-min)/n	#interval width

#function used to map a value to the intervals
#hist(x,width)=width*floor(x/width)+width/2.0

#set terminal postscript size 8.00,10.00 enhanced color font 'Helvetica,10' 
#set terminal postscript enhanced color eps font 'Helvetica,15' 
set terminal pdf transparent enhanced color font 'Helvetica,9' 
set output outfile

#set terminal png
#set output outdir."/".tool.".png"
set xrange [min:max]
set yrange [0:*]

#to put an empty boundary around the
set xtics min,(max-min)/5,max
#set boxwidth width*1.0
#set boxwidth 0.9 absolute

#set style data histogram
#set style data linespoints
#set style histogram rowstacked
#set style histogram clustered gap 1
#set style fill solid 0.5 border -1
#set style fill solid 0.25 noborder
set style fill transparent solid 0.5 border -1

#binwidth=5
set boxwidth binwidth*1.0

#function used to map a value to the intervals
bin(x,width)=width*floor(x/binwidth) + binwidth/2.0

#set style line 1 lt 1 lw 3 pt 3
#set style fill solid 0.5	#fillstyle
set tics out nomirror
#set title "Exome sample K24510-88962 (known indels -dbSNP 135)"
set xlabel "VAF (%)"
set ylabel "Frequency"
set xtic 20

#set multiplot layout 5, 1 title ""
#set multiplot layout 1, 3 title "Indel size distribution"
set tmargin 2

#set key box
#set key at 100, 10000
show key

set title "Homopolyers VAF distribution"
#set label 1 "$tool " at -80,400 center

#plot infile u (bin(100*($2/($1+$2)),binwidth)):(1.0) smooth freq w boxes title ""
plot infileA u (bin(100*($2/($1+$2)),binwidth)):(1.0) smooth frequency with boxes fs transparent solid  title "A", \
infileT u (bin(100*($2/($1+$2)),binwidth)):(1.0) smooth frequency with boxes fs transparent solid  title "T", \
infileG u (bin(100*($2/($1+$2)),binwidth)):(1.0) smooth frequency with boxes fs transparent solid  title "G", \
infileC u (bin(100*($2/($1+$2)),binwidth)):(1.0) smooth frequency with boxes fs transparent solid  title "C"
#unset multiplot
