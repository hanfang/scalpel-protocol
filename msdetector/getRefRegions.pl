#!/usr/bin/perl

#
# Script to generate fasta file input for quickscan
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

#parameters:
my $help;
my $VERBOSE=0;
my $REF = "";
my $VCF = "";
my $RADIUS = 50;

#variables:
my %genome;
my $start_time = time;
my $argcnt = scalar(@ARGV);

#####################################################
#
# Message about this program and how to use it
#
sub usage() {
	
print STDERR <<END;

usage: getRefRegions.pl [OPTIONS]

    --help    : this (help) message
    --verbose : verbose mode

    --vcf     : VCF file of mutations
    --ref     : reference file in FASTA format
    --radius  : radius of extension

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
	    
	# required parameters
    'vcf=s'     => \$VCF,
    'ref=s'     => \$REF,
	
	#optional 
    'radius=i'  => \$RADIUS,
    
) or usage();

# load the genome in fasta format
##########################################
sub loadGenomeFasta {
    
    print STDERR "Loading genome from FASTA file...";
    
    my $file = $_[0];
    my $genome = $_[1];
    
    open FASTAFILE, "< $file" or die "Can't open $file ($!)\n";

    my $header = "";
    my $SDNA = "";
    my $chr = "";
    my $cnt = 0;
	while ( <FASTAFILE> ) {
		chomp;
		#print "$_\n";
		if($_ =~ m/>/) {
			$cnt++;
			$header = $_;
			#print "$header\n";
			if($chr ne "") {
				#print "$chr\n";
				$genome->{$chr}->{seq} = $SDNA;
			}
            
			$header =~ s/^>//; # remove ">"
			$header =~ s/\s+$//; # remove trailing whitespace

			my ($label, $tmp) = split / /, $header, 2;
			$chr = $label;
			#if ($label =~ /^chr/) { $chr = substr($label,3); } # update chromosome label
			$SDNA = "";# clear out old sequence
		}
		else { 
			s/\s+//g; # remove whitespace
			$SDNA .= $_; 
		}
	}
    close FASTAFILE;

    if($chr ne "") { # handle last sequence
        #print "$chr\n";
        $genome->{$chr}->{seq} = $SDNA;
    }  

    #if($VERBOSE) {
    print STDERR "$cnt sequences.\n";
    #}
}


# get sequences from reference and 
# write them to output
##########################################
sub getSequences {
	
    print STDERR "Extract sequences...";
    
    my $file = $_[0];
	my $radius = $_[1];
    
    open VCFFILE, "< $file" or die "Can't open $file ($!)\n";
	
	my $cnt = 0;
    while ( <VCFFILE> ) {
        chomp;
        #print "$_\n";
        next if($_ =~ /^#/); # skip comments
   		$cnt++;
		
        my @records = split "\t", $_;
		my $chr = $records[0];
		# shift right by 1 bp to adjust for the VCF coordinate system 
		my $left = $records[1] - $radius + 1;
		my $right = $records[1] + $radius + 1;
		
		print ">$chr:$left-$right\t$_\n";
		die "Undefined sequence ($chr)\n" if (!exists($genome{$chr}));
		my $seq = substr($genome{$chr}->{seq}, $left-1, $right-$left+1);

		for(my $i=0; $i<length($seq); $i+=80) {
			my $str = substr($seq,$i,80);
			print "$str\n";
		}
	}
    print STDERR "..$cnt sequences.\n";
}

# print total time elapsed 
##########################################
sub elapsedTime {
	my $time_taken = $_[0];
	my $tool = $_[1];

	my $hours = $time_taken / 3600;
	my $days = $hours / 24; 
	$hours = $hours % 24;
	my $seconds = $time_taken % 3600;
	my $minutes = $seconds / 60;
	$seconds    = $seconds % 60;

	$days = floor($days);
	$hours = floor($hours);
	$minutes = floor($minutes);
	$seconds = floor($seconds);

	print STDERR "$tool elapsed time: $days day(s) $hours hour(s) $minutes minute(s) $seconds second(s)\n";
}

#####################################################
#
## do the job
##########################################

my $datestring = localtime();
if($VERBOSE) {
	print STDERR "Local date and time: $datestring\n";
}

#print STDERR "argcnt = $argcnt\n";
usage() if ($argcnt < 1);
usage() if ($help && ($argcnt==1));	
	
if( ($VCF eq "") && ($REF eq "")) { 
	print STDERR "Command required!\n";
	usage(); 
}
else {

	#load genome from fasta file
	loadGenomeFasta($REF, \%genome);
	getSequences($VCF, $RADIUS);
}

##########################################

my $time_taken = time - $start_time;

if($VERBOSE) {
	elapsedTime($time_taken, "Total");
}