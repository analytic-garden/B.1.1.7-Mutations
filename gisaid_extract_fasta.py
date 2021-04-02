#!/home/bill/anaconda3/envs/bio/bin/python
# -*- coding: utf-8 -*-
"""
ract_fasta.py
    Save sequences from MSA that are from a specified
    country of exposure or lineage.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_02_07

requires:
    Fasta file headers must have been reformatted 
    with gisaid_reformat_fasta_headers.py.
"""

import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def GetArgs():
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)
                
        parser = Parser(description='Save sequences from fasta file that are from a specified country or lineage.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input file (required). Fasta file headers must have been reformatted with gisaid_reformat_fasta_headers.py',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = False,
                            help = 'Output fasta file. Default: write to screen',
                            type = str)
        parser.add_argument('-c', '--country',
                            required = False,
                            help = 'Country of exposure',
                            type = str)
        parser.add_argument('-l', '--lineage',
                            required = False,
                            help = 'Pangolin lineage',
                            type = str)
        parser.add_argument('-r', '--replace',
                            required = False,
                            default = False,
                            help = "Replace U's with T's in sequence",
                            action = 'store_true')

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def write_record(rec, file, replace_u = False):
    """
    Write a single SeqRecord to an open file handle
    Optionally, convery U's in sequence to T

    Parameters
    ----------
    rec : BioPython SeqRecord
        A SeqRecord read by SeqIO.
    file : file handle
        Handle of a file open for writing.
    replace_u : bool, optional
        Replace U's in sequence by T. 
        The default is False.

    Returns
    -------
    None.

    """
    if not replace_u:
        SeqIO.write(rec, file, 'fasta')
    else:
        seq = str(rec.seq).replace('U', 'T')
        new_rec = SeqRecord(id = rec.id,
                            seq = Seq(seq),
                            description = rec.description)
        SeqIO.write(new_rec, file, 'fasta')
        
def main():
    args = GetArgs()
    msa_file = args.input_file
    country = args.country
    lineage = args.lineage
    out_file = args.output_file
    replace_u = args.replace
    
    if country is None and lineage is None:
        sys.exit()
       
    if not out_file is None:
        f = open(out_file, 'w')
    else:
        f = sys.stdout
        
    for rec in SeqIO.parse(msa_file, 'fasta'):
        items = rec.description.split('|')
        if country is not None and lineage is not None:
            if items[5] == country and items[4] == lineage:
                write_record(rec, f, replace_u)
        elif country is not None:            
            if items[5] == country:
                write_record(rec, f, replace_u)
        else:
            if items[4] == lineage:
                write_record(rec, f, replace_u)
            
    if not out_file is None:            
        f.close()

if __name__ == "__main__":
    main()
