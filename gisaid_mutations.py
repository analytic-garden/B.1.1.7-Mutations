#!/home/bill/anaconda3/envs/bio/bin python
# -*- coding: utf-8 -*-
"""
gisaid_mutations.py
    Find aligned columns with significant mutations 
    i.e. changes from reference seq.
    
@author: Bill Thompson
@license: GPL 3
@copyright: 2020_02_05

Note: This code contains a couple of big, ugly function that will 
need to be refactored.
"""

import sys
import argparse
from collections import Counter
from collections import defaultdict
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

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
                            help = 'Output CSV file (required).',
                            type = str)
        parser.add_argument('-g', '--refseq_gb_file',
                            required = False,
                            help = 'GenBank file for Wuhan reference sequence.',
                            type = str,
                            default = './NC_045512.gb')
        parser.add_argument('-q', '--min_column_quality',
                            required = False,
                            help = 'Minimum quality to be accepted. Quality is the percent of A, C, G, T in column.',
                            default = 0.9,
                            type = float)
        parser.add_argument('-c', '--consensus_cutoff',
                            required = False,
                            help = 'The upper value for variation in an alignment column.',
                            default = 0.99999,
                            type = float)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def read_alignment_file(filename):
    """
    read_alignment_file - read a collection of aligned sequences in FASTA format

    arguments:
        filename - the name of the FASTA file

    returns:
        a tuple:
            align_dict -a dictionary. The keys are sequence IDs. The values are BioPython SeqRecords.
            alignment - an AlignIO alignment object
    """
    alignment = AlignIO.read(open(filename), "fasta")

    align_dict = dict()
    for record in alignment:
        align_dict[record.id] = record

    return (align_dict, alignment)

def get_col_range(align, min_quality = 0.99):
    """

    Parameters
    ----------
    align : a BioPython Bio.Align.MultipleSeqAlignment object
    min_quality : a float
        the minimum column quality
        quality is the percent of A, C, G, T in column

    Returns
    -------
    start, end : ints
        the start and ending positions of positions with sufficient quality

    """    
    start = 0
    for col in range(0, align.get_alignment_length()):
        c = Counter(align[:, col]).most_common(4)
        qual = sum([x[1] for x in c if x[0] in ['A', 'C', 'G', 'T']])/len(align[:, col])
        if qual > min_quality:
            start = col
            break
        
    end = 0
    for col in range(align.get_alignment_length()-1, -1, -1):
        c = Counter(align[:, col]).most_common(4)
        qual = sum([x[1] for x in c if x[0] in ['A', 'C', 'G', 'T']])/len(align[:, col])
        if qual > min_quality:
            end = col
            break
    
    return start, end
    
def ref_pos_to_alignment(align, ref_seq):
    """
    Create a map of aligned positions to reference sequence positions.

    Parameters
    ----------
    align : dictionary
        An alignment dictioinary. Key is record ID.
    ref_seq : A BioPython Seq object
        A Seq object containg the reference sequence

    Returns
    -------
    pos_map : dictionary
        key = aligned column position, 
        value = reference coulum  position
    """
    pos_map = dict()
    ref_pos = 0    # 0 based
    for pos in range(len(ref_seq)):
        if ref_seq[pos] in ['A', 'C', 'G', 'T']:
            pos_map[pos] = ref_pos
            ref_pos += 1
        else:
            pos_map[pos] = -1

    return pos_map

def format_qual(qualifiers, id):
    """
    Format GenBank records qulaifiers.
    GenBank qualifiers contain commas. Since we are writing to a CVS
    we surround strings containing commas with quotes.
    See BioPython SeqIO docs for info about qualifiers.

    Parameters
    ----------
    qualifiers : dictionary
        key = GenBank ID, value =  GenBank qualifier
    id : str
        A GenBank ID

    Returns
    -------
    A string
        The qualifier string, possibly surrounded by quotes, 
        or an empty stringDESCRIPTION.
    """
    if not id in qualifiers:
        return ''

    qual = qualifiers[id][0]
    if ',' in qual:
        qual = '"' + qual + '"'

    return qual

def count_mutations(align, pos_map, ref_seq,
                    start=130, end=29840, consensus_cutoff=0.8):
    """
    Count the variation in an alignment
    This is a big ugly function that counts variation in columns
    from the reference. It returns a collection  of data for the 
    mutation CSV file.

    Parameters
    ----------
    align : A BioPython MultipleSeqAlignment object
        DESCRIPTION.
    pos_map : dictionary
        A map of aligned positions to reference positions generated by
        ref_pos_to_alignment()
    ref_seq : BioPython GenBank object
        Genbank recor of reference sequence.
    start : int, optional
        Start column in the alignment. The default is 130.
    end : int, optional
        Ending column in the alignment. The default is 29840.
        DESCRIPTION. The default is 29840.
    consensus_cutoff : float, optional
        The upper value for variation in an alignment column. 
        The default is 0.8.

    Returns
    -------
    dictionary
        alignment columns, reference columns - list of column numbers in alignment and reference sequence
        A,C,G,T,insertions - list the count of nucleutides and insertion in aligned columns
        other - list of the counts of everything else (N's ambiguous calls, etc.)                                                       consensus - a list consensus nucleotide for each column
        'consensus %' - a list consensus percent of consensus nucleotide for each column
        alt - a list of the next popular nucleotide in each column
        codons - list of codons containing consensus nucleotide
        aas - list amino acid of codon containing consensus nucleotide
        alt_codons - a list of the codons containing the next popular nucleotide 
        alt_aas - a list of the amino acids of the codons containing the next 
            popular nucleotide. 
        products, notes - protein and info GenBank listed by GenBank.
    """
    # helper function
    def remove_non_nucs(clist):
        """
        Remove non nucleotide characters from list of tuples 
        from Counter().

        Parameters
        ----------
        clist : list
            A list of tuples created by Counter(align[:, col])

        Returns
        -------
        clist : list
            A list of tuples with non-nucleotides removed.

        """
        remove_list = []
        for i in range(len(clist)):
            if clist[i][0] not in ['A', 'C', 'G', 'T', '-']:
                remove_list.append(i)
                
        for i in sorted(remove_list, reverse = True):
            clist.pop(i)
            
        return clist
        
    rows = align.__len__()

    align_cols = []
    ref_cols = []
    nts = {'A': [], 'C': [], 'G': [], 'T': [], '-': []}
    other = []
    consensus = []
    consensus_pct = []
    alt = []
    ref = []
    
    # run through columns and find coulmns with significant variation
    for col in range(start, end):          
        c = Counter(align[:, col])
        if len(c) > 1:
            common = c.most_common()
            
            if common[0][0] in ['N', '-']:
                continue   # ignore colums that are mostly N's or gaps
                
            # ignore non nucleotides
            common = remove_non_nucs(common)
            if len(common) <= 1:
                continue
            
            # mutations are variations from the reference
            if align[0, col] == common[0][0]:
                continue
            
            # if common[1][0] in ['N', '-']:  # ignore N's
            #     continue

            # count nucleotides
            counts = defaultdict(int)
            for nt, count in c.items():
                counts[nt] += count

            for nt in nts.keys():
                nts[nt].append(counts[nt])
            other.append(rows - sum([counts[k] for k in nts.keys()]))

            # get the consensus and next best
            consens = common[0][0]
            consensus.append(consens)
            try:
                consensus_pct.append((counts[consens] / sum([counts[k] for k in nts.keys()])) * 100)
            except:
                print('Invalid column:', col)
                sys.exit(1) 
                
            ref.append(align[0, col])
            alt.append(consens)
            
            align_cols.append(col+1)
            ref_cols.append(pos_map[col]+1)

    codons = []
    aas = []
    alt_codons= []
    alt_aas = []
    notes = []
    products = []
    refs = []
    aa_pos = []
    feat_type = []
    mutations = []
    new_nucs = []
    # find the codons and amino acids
    for position, pct, alt, consens in zip(ref_cols,
                                           consensus_pct,
                                           alt,
                                           ref):
        if pct/100 < consensus_cutoff: # only use significan columns
            seq_info = find_location_feature(ref_seq, position-1, alt, consens)
            feat_type.append(seq_info['feature_type'])
            
            if seq_info['feature_type'].find('UTR') == -1:
                aas.append(str(seq_info['aa']))
                codons.append(str(seq_info['codon']))
                alt_codons.append(str(seq_info['alt_codon']))
                alt_aas.append(str(seq_info['alt_aa']))
                mutations.append(str(seq_info['aa']) + \
                                 str(int(seq_info['aa_pos'])+1) + \
                                 str(seq_info['alt_aa']))
                new_nucs.append(str(seq_info['new_nuc']))
            else:
                aas.append('')
                codons.append('')
                alt_codons.append('')
                alt_aas.append('')
                mutations.append('')
                new_nucs.append(str(seq_info['new_nuc']))
                
            notes.append(format_qual(seq_info['qualifiers'], 'note'))
            products.append(format_qual(seq_info['qualifiers'], 'product'))
            refs.append(seq_info['ref_nuc'])
            aa_pos.append(seq_info['aa_pos']+1)
        else:
            codons.append('')
            aas.append('')
            alt_codons.append('')
            alt_aas.append('')
            notes.append('')
            products.append('')
            refs.append('')
            aa_pos.append('')
            feat_type.append('')
            mutations.append('')
            new_nucs.append('')

    return {'alignment columns': align_cols,
            'reference positions': ref_cols,
            'A': nts['A'],
            'C': nts['C'],
            'G': nts['G'],
            'T': nts['T'],
            'insertions': nts['-'],
            'other': other,
            'ref_nucleotide': refs,
            'consensus_nucleotide': consensus,
            'consensus %': consensus_pct,
            'alt_nucleotide': new_nucs,
            'feature_type': feat_type,
            'codon': codons,
            'ref_aa': aas,
            'alt_codons': alt_codons,
            'alt_aas': alt_aas,
            'aa_pos': aa_pos,
            'mutation': mutations,
            'product': products,
            'notes': notes
            }

def find_location_feature(ref_seq, position, alt, consens):
    """
    Get features at position from reference sequence.
    This is a helper function for count_mutations(). It should only
    by called by that function.

    Parameters
    ----------
    ref_seq : Bio.Seq GenBank record
        GenBank record for reference sequence.
    position : int
        Position in aligned sequence.
    alt : str
        An alternate nucleotide, most likely the most common
        nucleotide.
    consens : str
        Nucleotide from reference sequence.

    Returns
    -------
    dictionary
        conserved_position - the input position,
        start, end - start and end of feature
        feature_trype - a string, the feature type gene, UTR etc.
        codon_start, codon_end - atrt and end of codon
        codon - the codon,
        aa - the aa at position
        alt_codon, alt_aa - codon and amino acid if alt is substituted at position
        seq - sequence of feature
        qualifiers - feature qualifiers See Bio.Seq docs for details

    requires:
        position - must be in range of reference sequence length
        alt - A, C, G, or T
    """
    feat = ''
    for f in ref_seq.features:
        if position in f:
            feat = f
            
    feat_type = feat.type
       
    for p in feat.location.parts:
        start = p.start.real
        end = p.end.real
        if position >= start and position < end:
            break;
            
    sequence = ref_seq.seq[start:end]
        
    # start = feat.location.start.real
    # end = feat.location.end.real
    # sequence = ref_seq[start:(end+1)]
                       
    offset = position - start
    codon_pos = offset % 3
    codon_start = offset - codon_pos
    codon_end = codon_start + 2
    # codon = feat.extract(ref_seq).seq[codon_start:(codon_end+1)]
    codon = sequence[codon_start:(codon_end+1)]
    aa = codon.translate()
    
    ref_nuc = sequence[offset]
    new_nuc = alt
    if alt == ref_nuc:
        new_nuc = consens
    aa_pos = codon_start / 3
    
    seq_list = list(str(sequence))
    seq_list[offset] = new_nuc
    new_seq = SeqRecord(Seq(''.join(seq_list)))
    alt_codon = new_seq.seq[codon_start:(codon_end+1)]
    alt_aa = alt_codon.translate() if new_nuc in ['A', 'C', 'G', 'T'] else ''
    
    return {'conserved_position': position,
            'start': start,
            'end': end,
            'feature_type': feat_type,
            'codon_start': codon_start,
            'codon_end': codon_end,
            'codon': codon,
            'aa': aa,
            'alt_codon': alt_codon,
            'alt_aa': alt_aa,
            'ref_nuc': ref_nuc,
            'aa_pos': aa_pos,
            'seq': sequence,
            'qualifiers': feat.qualifiers,
            'new_nuc': new_nuc}

                    
def count_the_mutations(ref_seq_gb_file,
                        msa_file,
                        ref_id,
                        min_col_quality = 0.90,
                        consensus_cutoff = 0.99999):
    """
    Count mutations in msa_file.

    Parameters
    ----------
    ref_seq_gb_file : str
        The name of a GenBank file containing reference genome.
        The sequence ID is not the same as the GISAID ID, but
        the sequences are identical.
    msa_file : str
        The GISAID MSA file anme.
    ref_file : str
        A GISAID alignment containing a single FASTA record of
        teh aligned Wuhan reference sequence.
    min_col_quality : float, optional
        The minimum column quality. The default is 0.90.
        Quality is the percent of A, C, G, T in column.
    consensus_cutoff : float, optional
        The upper value for variation in an alignment column. 
        The default is 0.99999.

    Returns
    -------
    mutations : dictionary
        See count_mutations.count_mutations for a description.

    """
    ref_seq_gb = SeqIO.read(ref_seq_gb_file, 'genbank')
    
    align_dict, align = read_alignment_file(msa_file)
    ref_seq = align_dict[ref_id]
    
    start, end = get_col_range(align, min_col_quality)
    
    pos_map = ref_pos_to_alignment(align_dict, ref_seq.seq)    
    del align_dict  # remove this to save a bit of memory.
    
    mutations = count_mutations(align, 
                                pos_map,
                                ref_seq_gb,
                                start,
                                end, 
                                consensus_cutoff)
    
    return mutations                     
                    
def main():
    args = GetArgs()
    msa_file = args.input_file
    mutation_csv_file = args.output_file
    ref_seq_gb_file = args.refseq_gb_file
    min_col_quality = args.min_column_quality
    consensus_cutoff = args.consensus_cutoff
    
    ref_id = 'Wuhan/WIV04/2019'
        
    mutations = count_the_mutations(ref_seq_gb_file, 
                                    msa_file, 
                                    ref_id,
                                    min_col_quality = min_col_quality,
                                    consensus_cutoff = consensus_cutoff)
    
    df = pd.DataFrame(mutations)
    df = df.sort_values('reference positions')
    df.to_csv(mutation_csv_file, index=False)
    
if __name__ == "__main__":
    main()
