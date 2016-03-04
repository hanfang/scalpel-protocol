#!/usr/bin/env python

## ----------------------
# A python module for filtering Scalpel single sample VCF files
## ----------------------
# version 0.2.0 2015-09-17
# Han Fang (hanfang.cshl@gmail.com)
# Cold Spring Harbor Laboratory
## ----------------------
from __future__ import print_function
import os
import sys
import argparse
from io import open
import vcf

## the class for filtering scalpel vcf file 
class filter_vcf:
    'filter based on parental coverage, chi2 score, and alternative allele coverage'
    def __init__(self, in_vcf, mc, ac, zyg, cr, chi2, fp, out_vcf):
        self.in_vcf   = in_vcf
        self.mc = mc
        self.ac = ac
        self.zyg = zyg
        self.cr = cr
        self.chi2 = chi2
        self.fp = fp
        self.out_vcf = out_vcf
        self.sub_vcf = list()
    
    ## call function for reading and filtering 
    def __call__(self):
        self.vcf_reader = vcf.Reader(open(self.in_vcf, 'r'))
        
        if self.zyg == "both":
            for record in self.vcf_reader:
                if (record.INFO['MINCOV'] >= self.mc) & (record.INFO['ALTCOV'] >= self.ac) & (record.INFO['COVRATIO'] >= self.cr) & (record.INFO['CHI2'] <= self.chi2) & (record.INFO['FISHERPHREDSCORE'] >= self.fp):
                    self.sub_vcf.append(record)        
            return(self.sub_vcf)

        else:
            for record in self.vcf_reader:
                if (record.INFO['ZYG'] == self.zyg) & (record.INFO['MINCOV'] >= self.mc) & (record.INFO['ALTCOV'] >= self.ac) & (record.INFO['COVRATIO'] >= self.cr) & (record.INFO['CHI2'] <= self.chi2) & (record.INFO['FISHERPHREDSCORE'] >= self.fp):
                    self.sub_vcf.append(record)
            return(self.sub_vcf)

    ## write function for exporting vcf files
    def __write__(self):
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.out_vcf, 'w'), self.vcf_reader)
        elif int(sys.version_info.major) == 2 :
            vcf_writer = vcf.Writer(open(self.out_vcf, 'wb'), self.vcf_reader)

        for record in self.sub_vcf:
            vcf_writer.write_record(record)

## the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='single-vcf-filter.py - a script for filtering Scalpel single sample vcf')
    parser.add_argument("-i", help="input vcf file [REQUIRED]",required=True)
    parser.add_argument("-mc", type=int, default=0, help="keep varaints if MINCOV (minimum Kmer coverage for the indel) is greater than [INT], Default: 0")
    parser.add_argument("-ac", type=int, default=0, help="keep varaints if ALTCOV (Kmer coverage for the ref allele) is greater than [INT], Default: 0")
    parser.add_argument("-zyg", type=str, default="both", help="keep variants that ZYG is equal to [het/hom], or keep both types if not specified, Default: both")
    parser.add_argument("-cr",  type=float, default=0.0, help="keep varaints if COVRATIO (coverage ratio for the indel) is greater than [FLOAT], Default: 0.0")
    parser.add_argument("-chi", type=float, default=100, help="keep varaints if CHI2 square score is smaller than [float], Default: 100")
    parser.add_argument("-fp",  type=float, default=0, help="keep varaints if FISHERPHREDSCORE is greater than [FLOAT], Default: 0")
    parser.add_argument("-o",   help="output vcf file [REQUIRED]",required=True)

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        print ("[help] please use [-h] for details")
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.o!=None):
        foo = filter_vcf(args.i, int(args.mc), int(args.ac), str(args.zyg), float(args.cr), float(args.chi), float(args.fp), args.o )
        foo.__call__()
        foo.__write__()

    ## print usage message if any argument is missing
    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
