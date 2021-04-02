#!/home/bill/anaconda3/envs/bio2/bin/python
# -*- coding: utf-8 -*-
"""
gisaid_extract_ref_seq.py
        Get Wuhan reference sequence for MSA.
        Fasta file headers must have been reformatted with 
        gisaid_reformat_fasta_headers.py.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_03_23
"""
import argparse
import sys
from Bio import SeqIO

def GetArgs():
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Find aligned columns with significant mutations ')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input MSA file (required). Fasta file headers must have been reformatted with gisaid_reformat_fasta_headers.py',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = True,
                            help = 'Output Fasta file (required).',
                            type = str)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def main():
    args = GetArgs()
    input_file = args.input_file
    output_file = args.output_file
    
    ref_id = 'Wuhan/WIV04/2019'
    
    f = open(output_file, 'w')
    for rec in SeqIO.parse(input_file, 'fasta'):
        if rec.id == ref_id:
            SeqIO.write(rec, f, 'fasta')
            break
    f.close()
    
if __name__ == '__main__':
    main()

