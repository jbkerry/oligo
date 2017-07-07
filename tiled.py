'''The `tiled` module designs oligos for multiple adjacent
restriction fragments across a specified region of a chromosome, or for
the entire chromosome. If run in FISH mode (restriction fragment-
independent), it will generate end-to-end oligos based on a user-defined
step size.

'''

import re

from Bio import SeqIO

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}
    
def gen_oligos_capture(fa, chromosome, enzyme='DpnII', oligo=70, region=''):
    '''Generates fasta file containing the oligos for a specific
    restriction enzyme
    
    Parameters
    ----------
    fa : str
        Path to reference genome fasta
    chromosome : str
        Chromosome number/letter e.g. 7 or X
    enzyme : {'DpnII', 'NlaIII', 'HindIII'}, optional
        The enzyme for digestion, default=DpnII
    oligo : int, optional
        The length of the oligos to design (bp), default=70
    region : str
        The region of the chromosome to design oligos, must be in the
        format start-stop, e.g. 10000-20000; omit this option to design
        oligos over the entire chromosome
        
    Returns
    -------
    oligo_seqs : dict
        Contains all oligo sequences; key = oligo coordinate,
        value = oligo sequence

    '''
    
    chr_name = 'chr'+str(chromosome)
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    seq = seq_dict[chr_name].seq.upper()
    print('\t...complete\nGenerating oligos...')
        
    if not region:
        start = 0; stop = len(seq)
    else:
        start, stop = map(int, region.split('-'))
    
    p = re.compile(rs_dict[enzyme])
    pos_list = []
    for m in p.finditer(str(seq[start:stop])):
        pos_list.append(m.start()+start)
    
    cut_size = len(rs_dict[enzyme])
    
    oligo_seqs = {}
       
    for i in range(len(pos_list)-1):
        j = i + 1
        frag_len = pos_list[j]-pos_list[i]+cut_size  ## check this code
        if (frag_len>=oligo):
            l_start = pos_list[i]; l_stop = l_start+oligo
            r_stop = pos_list[j]+cut_size; r_start = r_stop-oligo
            
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

def gen_oligos_fish(fa, chromosome, step=70, oligo=70, region=''):
    '''Generates fasta file containing the oligos for a specific
    restriction enzyme
    
    Parameters
    ----------
    fa : str
        Path to reference genome fasta
    chromosome : str
        Chromosome number/letter e.g. 7 or X
    step : int, optional
        The step size, or how close adjacent oligos should be, default=70;
        for end-to-end oligos, set `step` to equal `oligo`
    oligo : int, optional
        The length of the oligos to design (bp), default=70
    region : str
        The region of the chromosome to design oligos, must be in the
        format start-stop, e.g. 10000-20000; omit this option to design
        oligos over the entire chromosome
        
    Returns
    -------
    oligo_seqs : dict
        Contains all oligo sequences; key = oligo coordinate,
        value = oligo sequence
    
    '''
    
    chr_name = 'chr'+str(chromosome)
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    seq = seq_dict[chr_name].seq.upper()
    print('\t...complete\nGenerating oligos...')
        
    if not region:
        start = 0; stop = len(seq)
    else:
        start, stop = map(int, region.split('-'))
    
    oligo_seqs = {}
    
    final_start = stop-oligo

    while start<=final_start:
        stop = start + oligo
        ol_seq = str(seq[start:stop])
        oligo_seqs['{}:{}-{}-000-000-X'.format(chr_name, start, stop)] = ol_seq
        start+=step  
    
    print('\t...complete')
    return oligo_seqs
    
def split_fa():
    f = open('oligo_seqs.fa')
    file_content = f.readlines()
    split = 40000
    start = 1
    stop = 20000
    for lines in range(0, len(file_content), split):
        output_data = file_content[lines:lines+split]
        with open('{}-{}.fa'.format(start, stop),'w') as f_out:
            f_out.write(''.join(output_data))
        start+=20000; stop+=20000
        
    print('Split files per 20,000 oligos')
    return

if __name__ == '__main__':
    
    import tools
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fish',
        action = 'store_true',
        help = 'Run in FISH mode (restriciton enzyme independent).',
        required = False,
    )
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
        help = 'Name of restriction enzyme, default=DpnII. Omit this option ' \
               'if running in FISH mode (--fish)',
        default = 'DpnII',
        required = False,
    )
    parser.add_argument(
        '-t',
        '--step_size',
        type = int,
        help = 'Step size (in bp) between adjacent oligos when running in ' \
               'FISH mode (--fish), default=70. Omit this option if you are ' \
               'not using the --fish flag',
        default = 70,
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
               'format \'start-stop\' e.g. \'10000-20000\'. Omit this ' \
               'option to design oligos across the entire chromosome.',
        required = False,
    )
    parser.add_argument(
        '-s',
        '--star_index',
        type = str,
        help = 'Path to STAR index directory. Omit this option if running ' \
               'with BLAT (--blat)',
        required = False,
    )
    parser.add_argument(
        '--blat',
        action = 'store_true',
        help = 'Detect off-targets using BLAT instead of STAR.',
        required = False,
    )
    parser.add_argument(
        '--test_fasta',
        action = 'store_true',
        help = argparse.SUPPRESS,
        required = False,
    )
    
    args = parser.parse_args()
    
    if not args.blat and not args.star_index:
        msg = '-s/--star_index argument is required if --blat is not selected'
        parser.error(msg)
    
    if args.fish:
        pass_seqs = gen_oligos_fish(
            fa = args.fasta,
            chromosome = args.chr,
            step = args.step_size,
            oligo = args.oligo,
            region = args.region,
        )
    else:
        pass_seqs = gen_oligos_capture(
            fa = args.fasta,
            chromosome = args.chr,
            enzyme = args.enzyme,
            oligo = args.oligo,
            region = args.region,
        )
    tools.write_oligos(oligo_seqs=pass_seqs)
    if not args.test_fasta:
        tools.check_off_target(
            genome = args.genome,
            fa = args.fasta,
            s_idx = args.star_index,
            blat=args.blat,
        )
        tools.get_density(blat=args.blat)
