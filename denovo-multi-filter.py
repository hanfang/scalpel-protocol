#!/usr/bin/env python

## ----------------------
# A python module for filtering Scalpel multisample VCF files
## ----------------------
# version 0.2.0 2015-09-18
# Han Fang (hanfang.cshl@gmail.com)
# Cold Spring Harbor Laboratory
## ----------------------
from __future__ import print_function
import os
import sys
import argparse
from io import open
import vcf

## the class for filtering scalpel multisample vcf file 
class filter_vcf:
    'filter based on parental coverage, chi2 score, and alternative allele coverage'
    def __init__(self, in_vcf, father, mother, affected, unaffected, aac, chi2, pc, out_vcf):
        self.in_vcf   = in_vcf
        self.father   = father
        self.mother   = mother
        self.affected = affected
        self.unaffected = unaffected
        self.aac = aac
        self.chi2 = chi2
        self.pc = pc
        self.out_vcf = out_vcf
        self.sub_vcf = list()
    
    ## call function for reading and filtering 
    def __call__(self):
        self.vcf_reader = vcf.Reader(open(self.in_vcf, 'r'))
        
        for record in self.vcf_reader:
            if ( (record.genotype(self.father)['DP'] >= self.pc) & (record.genotype(self.mother)['DP'] >= self.pc) & (record.genotype(self.affected)['AD'][1] >= self.aac) & (record.INFO['CHI2'] <= self.chi2) ):
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
    parser = argparse.ArgumentParser(description='denovo-multi-filter.py - a script for filtering Scalpel multi-sample vcf')
    parser.add_argument("-i", help="input vcf file [REQUIRED]",required=True)
    parser.add_argument("-f", help="sample id of the father [REQUIRED]",required=True)
    parser.add_argument("-m", help="sample id of the mother [REQUIRED]",required=True)
    parser.add_argument("-a", help="sample id of the affected child [REQUIRED]",required=True)
    parser.add_argument("-u", help="sample id of the unaffected child [REQUIRED]",required=True)
    parser.add_argument("-aac", type=int  , default=0, help="keep varaints if MINCOV (alternative allele coverage) is greater than [INT], Default: 0")
    parser.add_argument("-chi", type=float, default=100, help="keep varaints if CHI2 Square score is smaller than [float], Default: 100")
    parser.add_argument("-pc",  type=int  , default=0, help="keep varaints if the coverage in either parent is greater than [INT], Default: 0")
    parser.add_argument("-o",   help="output vcf file [REQUIRED]",required=True)

    ## check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        print ("[help] please use [-h] for details")
        sys.exit(1)
    else:
        args = parser.parse_args()

    ## process the file if the input files exist
    if (args.i!=None) & (args.f!=None) & (args.m!=None) & (args.a!=None) & (args.u!=None) & (args.o!=None):
        foo = filter_vcf(args.i, args.f, args.m, args.a, args.u, int(args.aac), float(args.chi), int(args.pc), args.o )
        foo.__call__()
        foo.__write__()

    ## print usage message if any argument is missing
    else:
        print ("[error]\tmissing argument")
        parser.print_usage()
