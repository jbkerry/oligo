#!/usr/bin/env python

from Bio import SeqIO

def gen_oligos(fa, bed, enzyme='DpnII', oligo=70):
    '''Generates oligos flanking restriction fragments that encompass the
    coordinates supplied in the bed file
    
    Parameters
    ----------
    fa: path to reference genome fasta
    bed: tab-delimited bed file containing a list of coordinates for viewpoints
        in the capture experiment.
        Must be in the format 'chr'\t'start'\t'stop'\t'name'
    enzyme: DpnII (GATC), NlaIII (CATG) or HindIII (AAGCTT), default=DpnII
    oligo: the length of the oligo to design (bp), default=70
    
    Returns
    -------
    
    '''
    
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    with open(bed) as w:
        for x in w:
            chr_name, start, stop, name = x.rstrip('\n').split('\t')
            
        
