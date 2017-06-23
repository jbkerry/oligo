#!/usr/bin/env python

import math
from Bio import SeqIO

def gen_oligos(fa, bed, oligo=70, step=10, max_dist=200):
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
    oligo_seqs: dictionary containing oligo coordinates as keys and oligo
        sequences as items
    
    '''
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    print('\t...complete\nGenerating oligos...')

    oligo_num = math.floor((max_dist-oligo)/step)
    oligo_seqs = {}
    
    with open(bed) as w:
        for x in w:
            chr_name, start, stop, name = x.rstrip('\n').split('\t')
            
            start, stop = map(int, (start, stop))
            seq = seq_dict[chr_name].seq.upper()
            
            l_stop = start-10
            l_start = l_stop-oligo
            r_start = stop+10
            r_stop = r_start+oligo
            
            counter=1
            while counter<=oligo_num:
                
                l_seq = seq[l_start:l_stop]
                r_seq = seq[r_start:r_stop]
                oligo_seqs['{}:{}-{}-000-000-X'.format(chr_name,
                                                       l_start,
                                                       l_stop)] = str(l_seq)
                oligo_seqs['{}:{}-{}-000-000-X'.format(chr_name,
                                                       r_start,
                                                       r_stop)] = str(r_seq)
                l_start-=step
                l_stop-=step
                r_start+=step
                r_stop+=step
                counter+=1
    
    print('\t...complete')
    return oligo_seqs