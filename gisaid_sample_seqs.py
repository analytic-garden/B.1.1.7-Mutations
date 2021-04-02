#!/home/bill/anaconda3/envs/bio2/bin/python
# -*- coding: utf-8 -*-
"""
gisaid_sample_seqs.py
        Sample a collectionf sequences fron GISAID sequence file.
        Fasta file headers must have been reformatted with 
        gisaid_reformat_fasta_headers.py.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_03_31
"""
import argparse
import sys
import random
from Bio import SeqIO

def GetArgs():
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Sample sequences form A GISAID FASTA file.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input FASTA file (required). Fasta file headers must have been reformatted with gisaid_reformat_fasta_headers.py',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = True,
                            help = 'Output Fasta file (required).',
                            type = str)
        parser.add_argument('-n', '--num_samples',
                            required = False,
                            help = 'Number of sequences to sample. (default = 1000)',
                            default = 1000,
                            type = int)
        
        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def main():
    args = GetArgs()
    input_file = args.input_file
    output_file = args.output_file
    num_samples = args.num_samples
    
    seqs = dict()
    for seq in SeqIO.parse(input_file, 'fasta'):
        seqs[seq.id] = seq
        
    seq_ids = random.sample(list(seqs.keys()), k=num_samples)
    
    f = open(output_file, 'w')
    for k in seq_ids:
        SeqIO.write(seqs[k], f, 'fasta')
    f.close()
    
if __name__ == '__main__':
    main()