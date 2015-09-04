#!/bin/sh

#
# Detects mutations which are within a STR (as defined below)
# Adds 3 columns to the VCF file with the status of the mutations
#
#
# version 1.0 2015-Feb-23
# Giuseppe Narzisi (gnarzisi@nygenome.org)
# New York Genome Center
#

# By default, it finds all microsatellites that are at least 8bp long (total length), 
# where the repeat sequence is between 1bp and 4bp, and is repeated at least 3 times.

#$ (echo ">A"; echo "CAAAAAAAAAAAAAAAAAAAAAAAT") > a.fa
#$ ./quickscan -f a.fa
#Processing sequences in a.fa...
#A len=25
#2 24 23 A 22+1 cAAAAAAAAAAAAAAAAAAAAAAAt 0
#finished in 8.9e-05s

SCRIPT=$(readlink -f "$0")
BINDIR=$(dirname "$SCRIPT")
REF=/data/NYGC/Resources/Indexes/bwa/human_g1k_v37.fasta

quickscan=$BINDIR/quickscan
getregions=$BINDIR/getRefRegions.pl
parse=$BINDIR/parseQiuickscan.pl


varFileList=""
radius=50
delta=1
units=3 # minimum number of units to report
strbp=8 # minimum length of tandem in bp 
maxulen=4 # max unit length 
minulen=1 # min unit length 
flkbp=10 # flanking bp to report


USAGE="Usage: `basename $0` [-hv] \n
\t -i <vcf>    VCF file with list of variants \n
\t -r <radius> sequence radius around the location to examine \n
\t -d <delta>  distance (bp) of STR from locus \n
\t -g <genome> genome file in FASTA format \n 
\t -u <units>  minimum number of units to report (default: $units) \n
\t -l <bp>     minimum length of tandem in bp (default: $strbp) \n
\t -x <len>    max unit length (default: $maxulen) \n
\t -m <len>    min unit length (default: $minulen) \n
\t -k <bp>     flanking bp to report (default: $flkbp) \n"

usage()
{
    echo -e $USAGE >&2
    exit 1
}

[ "$#" -lt 1 ] && usage

# Parse command line options.
while getopts hvi:r:d:g: OPT; do
	case "$OPT" in
		h)
			echo $USAGE
			exit 0
			;;
		v)
			echo "`basename $0` version 0.1"
			exit 0
			;;
		i)
			varFileList=$OPTARG
			;;
		r)
			radius=$OPTARG
			;;
		d)
			delta=$OPTARG
			;;
		g)
			REF=$OPTARG
			;;
		u)
			units=$OPTARG
			;;
		l)
			strbp=$OPTARG
			;;
		x)
			maxulen=$OPTARG
			;;
		m)
			minulen=$OPTARG
			;;
		k)
			flkbp=$OPTARG
			;;
		\?)
			# getopts issues an error message
			echo $USAGE >&2
			exit 1
			;;
	esac
done

# Remove the options we parsed above.
shift `expr $OPTIND - 1`

# We want at least one non-option argument.
# Remove this block if you don't need it.
#if [ $# -eq 0 ]; then
#    echo $USAGE >&2
#    exit 1
#fi

if [ ! -f $REF ]; then
    echo "Reference genome cannot be found" >&2
    exit 1
fi

if [ ! -f $quickscan.cc ]; then
    echo "Script: $quickscan.cc cannot be found" >&2
    exit 1
fi

if [ ! -f $getregions ]; then
    echo "Script: $getregions cannot be found" >&2
    exit 1
fi

if [ ! -f $parse ]; then
    echo "Script: $parse cannot be found" >&2
    exit 1
fi


# grab header
awk '$0~/^#/ {print}' $varFileList

# create list of regions to examine in fasta format
FASTA="${varFileList}.regions.fa"
$getregions --vcf $varFileList --ref $REF --radius $radius > $FASTA

# run microsatellites detector (quickscan) on regions file
# By default, it finds all microsatellites that are at least 8bp 
# long (total length), where the repeat sequence is between 
# 1bp and 4bp, and is repeated at least 3 times. 

QSOUT="${varFileList}.quickscan.out"
$quickscan -f $FASTA -u $units -l $strbp -x $maxulen -m $minulen -k $flkbp > $QSOUT
rm $FASTA

# parse regions file from quickscan
#MSOUT="regions.out"
$parse -i $QSOUT -r $radius -d $delta
rm $QSOUT
