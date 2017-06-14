#!/usr/bin/env python

import sys
import re
import argparse
from Bio import SeqIO
import tools

spe_dict = {'mm9': 'mouse',
            'mm10': 'mouse',
            'hg18': 'human',
            'hg19': 'human',
            'hg38': 'human'}

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}
    
def gen_oligos_capture(fa, chromosome, enzyme='DpnII', oligo=70, region=''):
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
    
    chr_name = 'chr'+str(chromosome)
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    seq = seq_dict[chr_name].seq.upper()
    print('\t...complete\nGenerating oligos...')
        
    if not region:
        start = 0; stop = len(seq)
    else:
        start, stop = list(map(int, region.split('-')))
    
    p = re.compile(rs_dict[enzyme])
    pos_list = []
    for m in p.finditer(str(seq[start:stop])):
        pos_list.append(m.start()+start)
    
    cut_size = len(rs_dict[enzyme])
    
    fa_w = open('oligo_seqs.fa', 'w')    
    for i in range(len(pos_list)-1):
        j = i + 1
        frag_len = pos_list[j]-pos_list[i]+cut_size
        if (frag_len>=oligo):
            l_start = pos_list[i]
            l_stop = l_start+oligo
            l_seq = seq[l_start:l_stop]
            
            r_stop = pos_list[j]+cut_size
            r_start = r_stop-oligo
            r_seq = seq[r_start:r_stop]
            
            fa_w.write('>{0}:{1}-{2}-{1}-{3}-L\n{4}\n'.format(chr_name,
                                                              l_start,
                                                              l_stop,
                                                              r_stop,
                                                              l_seq))
            if frag_len>oligo:
                fa_w.write('>{0}:{1}-{2}-{3}-{2}-R\n{4}\n'.format(chr_name,
                                                                  r_start,
                                                                  r_stop,
                                                                  l_start,
                                                                  r_seq))
    fa_w.close()
    
    print('\t...wrote oligos to oligo_seqs.fa')
    
def split_fa():
    f = open('oligo_seqs.fa')
    file_content = f.readlines()
    split = 40000
    start = 1
    stop = 20000
    for lines in range(0, len(file_content), split):
        output_data = file_content[lines:lines+split]
        with open('{0}-{1}.fa'.format(start, stop),'w') as f_out:
            f_out.write(''.join(output_data))
        start+=20000; stop+=20000
        
    print('Split files per 20,000 oligos')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f',
        '--fasta',
        type = str,
        help = 'Path to reference genome fasta.',
        required = True,
    )
    parser.add_argument(
        '-g',
        '--genome',
        type = str,
        help = 'Genome build e.g. \'mm10\' or \'hg38\'.',
        required = True,
    )
    parser.add_argument(
        '-c',
        '--chr',
        type = str,
        help = 'Chromosome number/letter on which to design the oligos.',
        required = True,
    )
    parser.add_argument(
        '-e',
        '--enzyme',
        type = str,
        help = 'Name of restriction enzyme, default=DpnII.',
        default = 'DpnII',
        required = False,
    )
    parser.add_argument(
        '-o',
        '--oligo',
        type = int,
        help = 'The size (in bp) of the oligo to design, default=70',
        default = 70,
        required = False,
    )
    parser.add_argument(
        '-r',
        '--region',
        type = str,
        help = 'The region in which to design the oligos; must be in the ' \
             'format \'start-stop\' e.g. \'10000-20000\'. Omit this option ' \
             'to design oligos across the entire chromosome.',
        required = False,
    )
    parser.add_argument(
        '-s',
        '--star_index',
        type = str,
        help = 'Path to STAR index directory, omit this if running with ' \
               'BLAT (--blat)',
        required = False,
    )
    parser.add_argument(
        '--blat',
        action = 'store_true',
        help = 'Detect off-targets using BLAT instead of STAR.',
        required = False,
    )
    args = parser.parse_args()
    
    gen_oligos_capture(
        fa = args.fasta,
        chromosome = args.chr,
        enzyme = args.enzyme,
        oligo = args.oligo,
        region = args.region
    )
    tools.check_off_target(
        species = spe_dict[args.genome.lower()],
        fa = args.fasta,
        s_idx = args.star_index,
        blat=args.blat
    )
    tools.get_density(blat=args.blat)
