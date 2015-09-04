#!/usr/bin/perl

#
# Script to parse the output of the quyickscan program
#
#
# version 1.0 2015-Feb-23
# Giuseppe Narzisi (gnarzisi@nygenome.org)
# New York Genome Center
#

use warnings;
use strict;
use POSIX;
use Getopt::Long;

my $help;
my $VERBOSE = 0;

# required arguments
my $INFILE;
my $R;
my $D=1;

# optional arguments
#my $TYPE;


my $version_num="1.0";
my $argcnt = scalar(@ARGV);

#####################################################
#
# Message about this program and how to use it
#
sub header() {
print STDERR <<END;

Program: parseQiuickscan
Version: $version_num
Contact: Giuseppe Narzisi <gnarzisi\@nygenome.org>
END
}

#####################################################
#
# Message about this program and how to use it
#
sub usage() {
	
header();
	
print STDERR <<END;

usage: parseQiuickscan --in <quickscan fasta file> [OPTIONS]

COMMAND:

    --help    : this (help) message
    --verbose : verbose mode

    Required parameters:
    --in      : list of STRs signatures detected by quickscan (in fasta format)
    --radius  : region radius 

    Optional parameters:
    --delta   : number of bp to check around the locus

END
exit;
}

#####################################################
#
# Command line options processing
#
GetOptions(
	'help!'    => \$help,
	'verbose!' => \$VERBOSE,

	# required parameters (single mode)
    'in=s'     => \$INFILE,
	'radius=i' => \$R,

	# optional paramters
	'delta=i'  => \$D,
	
) or usage();

#####################################################
#
# Initilization 
#
sub init()
{
	#print STDERR "argcnt = $argcnt\n";
	usage() if ($argcnt < 1);	
	usage() if ($help && ($argcnt==1));	
	
	if(!defined $INFILE) { 
		print STDERR "Command required!\n";
		usage(); 
	}
	else {		
		# run analysis
		go();
	}
}

# load the genome in fasta format
##########################################
sub parseGenomeFasta {
    
    print STDERR "Parse quickscan FASTA file...\n";
    
    my $file = $_[0];

    
    open FASTAFILE, "< $file" or die "Can't open $file ($!)\n";

    my $header = "";
    my $label = "";
	my $cnt = 0;
	my $status = "no";
	my $motif = "";
	my $L = "";
    while ( <FASTAFILE> ) {
        chomp;
        #print "$_\n";
        if($_ =~ m/>/) {
			$cnt++;
            $header = $_;
            #print "$header\n";
			 
			if($label ne "") { print STDOUT "$label\t$L\t$motif\t$status\n"; }
            
			# reset status and moitif
			$status="no";
			$motif="";
			$L = "";
			
            $header =~ s/^>//; # remove ">"
            $header =~ s/\s+$//; # remove trailing whitespace

            my ($hd, $tmp) = split / /, $header, 2;
            my ($reg, $mut) = split /\t/, $hd, 2;
			$label = $mut;
        }
        else { 
            #s/\s+//g; # remove whitespace
			# start end total_len tandem_unit complete_units+partial left_flank+repeat+right_flank
			# 51      67      17      A       17+0    actccatctcAAAAAAAAAAAAAAAAAgatgttaaca   0
            my ($start, $end, $len, $unit, $complete, $seq, $num) = split /\t/, $_, 7;
			$motif = $unit;
			$L=$len;
			# include microsatellite adjacent to the mutation
			if ( ($R+1)>=($start-$D) && ($R+1)<=($end+$D)) { $status="yes"; } 
			#print STDERR "$label\t$start\t$end\t$len\t$motif\t$status\n";
			#if ( ($R+1)>=$start && ($R+1)<=$end) { $status="yes"; } 
        }
    }
    close FASTAFILE;

    if($label ne "") { # handle last sequence
        print STDOUT "$label\t$L\t$motif\t$status\n";
	}

	#if($VERBOSE) {
		print STDERR "$cnt regions.\n";
	#}
}

#####################################################
#
# load annovare file
#
sub go{
	parseGenomeFasta($INFILE);
	#analysis();
}


#####################################################
#
# do the job
#
init();




