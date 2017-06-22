#!/usr/bin/env python

from Bio import SeqIO
import numpy as np
import re

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

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
    
    cut_sites = {}
    oligo_seqs = {}
    cut_size = len(rs_dict[enzyme])
    p = re.compile(rs_dict[enzyme])
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    for x in seq_dict:
        seq = str(seq_dict[x].seq).upper()
        cut_sites[x] = [m.start() for m in p.finditer(seq)]
        cut_sites[x] = np.array(cut_sites[x])
    print('\t...complete\nGenerating oligos...')
    
    with open(bed) as w:
        for x in w:
            chr_name, start, stop, name = x.rstrip('\n').split('\t')
            
            start, stop = map(int, (start, stop))
            seq = seq_dict[chr_name].seq.upper()
            
            l_start = cut_sites[chr_name][cut_sites[chr_name]<=start][-1]
            r_stop = cut_sites[chr_name][cut_sites[chr_name]>=start][0] + cut_size
            frag_len = r_stop - l_start
            
            if frag_len>=oligo:
                l_stop = l_start + oligo
                r_start = r_stop - oligo
                
                l_tup = (l_start, l_stop, l_start, r_stop)
                r_tup = (r_start, r_stop, l_start, r_stop)
                l_seq = seq[l_start:l_stop]
                r_seq = seq[r_start:r_stop]
                
                oligo_seqs['{}:{}-L'.format(chr_name,
                                           '-'.join(map(str, l_tup)))
                            ] = str(l_seq)
                if frag_len>oligo:
                    oligo_seqs['{}:{}-R'.format(chr_name,
                                           '-'.join(map(str, r_tup)))
                            ] = str(r_seq)
    
    print('\t...complete')
    return oligo_seqs
