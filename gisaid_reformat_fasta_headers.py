#!/home/bill/anaconda3/envs/bio/bin/python
# -*- coding: utf-8 -*-
"""
gisaid_reformat_fasta_headers.py
    Reformat fasta headers to make gisaid_epi_isl record id.
    New fasta header contains country of exposure, lineaage, and 
    gisaid_epi_isl as elements 4 through 6 of description.

@author: Bill Thompson
@license: GPL 3
@copyright: 2021_02_07
"""

import sys
import argparse
import pandas as pd
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
                
        parser = Parser(description='Reformat FASTA headers to include lineage, isl id, and country of exposure.')
        parser.add_argument('-i', '--input_file',
                            required = True,
                            help = 'Input file (required). Fasta file from GISAID.',
                            type = str)
        parser.add_argument('-m', '--meta_file',
                            required = True,
                            help = 'Meta file (required). GISAID meta file corresponding to FASTA file.',
                            type = str)
        parser.add_argument('-o', '--output_file',
                            required = False,
                            help = 'Output fasta file. Default: write to screen.',
                            type = str)
        parser.add_argument('-r', '--replace',
                            required = False,
                            default = False,
                            help = "Replace U's with T's in sequence.",
                            action = 'store_true')
        parser.add_argument('-n', '--min_Ns',
                            required = False,
                            default = 20,
                            help = "Minimum allowed number of N's in a row. (default = 20)",
                            type = int)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def main():
    args = GetArgs()
    msa_file = args.input_file
    meta_file = args.meta_file
    out_file = args.output_file
    replace_u = args.replace
    min_Ns = args.min_Ns
    
    ns = 'N' * min_Ns
    
    if out_file is not None:
        f = open(out_file, 'w')
    else:
        f = sys.stdout
    
    df = pd.read_csv(meta_file, 
                     sep = '\t', 
                     dtype={'location': str})
       
    for rec in SeqIO.parse(msa_file, 'fasta'):
        if str(rec.seq).find(ns) != -1:
            print("Too many N's:", rec.description, file = sys.stderr)
            continue
            
        items1 = rec.description.split('|')
        items2 = items1[0].split('/')
        new_id = '/'.join(items2[1:]).replace(' ', '')
        new_id = new_id.replace("'", '-')
        new_id = new_id.replace(',','')
        
        df2 = df[df['strain'] == new_id]
        
        if df2.size == 0:
            print('Missing ID in metafile:', rec.description, 
                  file = sys.stderr)
            continue
                
        new_desc = '|'.join([rec.description,
                             str(df2['pango_lineage'].values[0]),
                             str(df2['country_exposure'].values[0]),
                             str(df2['gisaid_epi_isl'].values[0])])
        
        if not replace_u:
            seq = rec.seq
        else:
            seq = Seq(str(rec.seq).replace('U', 'T'))

        new_rec = SeqRecord(seq = seq,
                            id = new_id,
                            description = new_desc)
        SeqIO.write(new_rec, f, 'fasta')
        
    if out_file is not None:            
        f.close()
                 
if __name__ == "__main__":
    main()
