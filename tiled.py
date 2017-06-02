#!/usr/bin/env python

import re
from Bio import SeqIO

org_dict = {'mm9': 'Mus_musculus',
            'mm10': 'Mus_musculus',
            'hg18': 'Homo_sapiens',
            'hg19': 'Homo_sapiens',
            'hg38': 'Homo_sapiens'}

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

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
            
        Output
        ------
        oligo_seqs.fa: a FASTA file containing sequences of all the oligos
        
        '''
        
        sequence = SeqIO.parse('/databank/igenomes/{0}/UCSC/{1}/Sequence/' \
            'WholeGenomeFasta/genome.fa'.format(org_dict[self.genome],
                                                self.genome), 'fasta')
        seq_dict = {}
        for seq_record in sequence:
            seq_dict[seq_record.name] = seq_record.seq.upper()
            
        if not region:
            start = 0
            stop = len(seq_dict[chromosome])
        else:
            start, stop = list(map(int, region.split('-')))
            
        p = re.compile(rs_dict[enzyme])
        pos_list = []
        for m in p.finditer(seq_dict[chromosome][start:stop]):
            pos_list.append(m.start()+start)
        
        fa = open('oligo_seqs.fa', 'w')    
        for i in range(len(pos_list)):
            j = i + 1
            frag_len = pos_list[j]-pos_list[i]
            if (j<len(pos_list)) & (frag_len>=oligo):
                l_start = pos_list[i]
                l_stop = l_start+oligo
                l_seq = seq_dict[chromosome][l_start:l_stop]
                
                r_stop = pos_list[j]
                r_start = r_stop-oligo
                r_seq = seq_dict[chromosome][r_start:r_stop]
                
                fa.write('>{0}:{1}-{2}-{1}-{3}-L\n{4}\n'.format(chromosome,
                                                                   l_start,
                                                                   l_stop,
                                                                   pos_list[j],
                                                                   l_seq))
                if fragment_length>oligo:
                    fa.write('>{0}:{1}-{2}-{3}-{2}-L\n{4}\n'.format(chromosome,
                                                                    r_start,
                                                                    r_stop,
                                                                    pos_list[i],
                                                                    r_seq))
        fa.close()
        
        return "Wrote oligos to oligo_seqs.fa"
