#!/usr/bin/env python

import re
from Bio import SeqIO

org_dict = {'mm9': 'Mus_musculus',
            'mm10': 'Mus_musculus',
            'hg18': 'Homo_sapiens',
            'hg19': 'Homo_sapiens',
            'hg38': 'Homo_sapiens'}

class Capture(object):
    '''Designs oligos for capture from adjacent restriction sites within a
    user-specified region
    
    Parameters
    ----------
    genome: hg18, hg19, hg38, mm9 or mm10
    blat: boolean, check off-targets using BLAT instead of STAR, default=False
    '''
    
    def __init__(self, genome, blat=False):
        self.genome = genome.lower()
        self.blat = blat
    
    def generate_oligos(self, chromosome, enzyme='DpnII', oligo=70, region=''):
        '''Generates fasta file containing the oligos for a specific
        restriction enzyme
        
        Parameters
        ----------
        chromosome: chromosome number/letter e.g. 7 or X
        enzyme: DpnII (GATC), NlaIII (CATG) or HindIII (AAGCTT), default=DpnII
        oligo: the length of the oligo to design (bp), default=70
        region: the region of the chromosome to design oligos, must be in the
            format start-stop, e.g. 10000-20000. Omit this option to design
            oligos over the entire chromosome
        
        '''
        
        _sequence = SeqIO.parse('/databank/igenomes/{0}/UCSC/{1}/Sequence/' \
                                'WholeGenomeFasta/genome.fa'.format(org_dict[self.genome], self.genome),
                                'fasta')
        _sequence_dict = {}
        for seq_record in _sequence:
            _sequence_dict[seq_record.name] = seq_record.seq.upper()
            
        if not region:
            _start_seq = 0
            _stop_seq = len(self._sequence_dict[chromosome])
        else:
            _start_seq, _stop_seq = region.split("-")
            try:
                _start_seq = int(_start_seq)
                _stop_seq = int(_stop_seq)
            except ValueError:
                print("Chromosome coordinates must be integers")
                raise
        
