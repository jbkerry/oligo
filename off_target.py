'''The `off_target` module designs oligos adjacent to user-supplied
coordinates for potential CRISPR off-target cleavage sites, to allow for
efficient pull-down of the edited region

'''

import math
import pickle
import sys

from Bio import SeqIO

def gen_oligos(fa, bed, oligo=70, step=10, max_dist=200):
    r'''Generates oligos adjacent to user-supplied coordinates for
    potential CRISPR off-target cleavage sites
    
    Parameters
    ----------
    fa : str
        Path to reference genome fasta
    bed : str
        Path to tab-delimited bed file containing a list of coordinates
        for viewpoints in the capture experiment; must be in the format
        'chr'\\t'start'\\t'stop'\\t'name'\\n
    oligo : int, optional
        The length of the oligo to design (bp), default=70
    step : int, optional
        The step size, or how close adjacent oligos should be,
        default=10; for end-to-end oligos, set `step` to equal `oligo`
    max_dist : int, optional
        The maximum distance away from the off-target site to design
        oligos to, default=200
    
    Returns
    -------
    oligo_seqs : dict
        Contains all oligo sequences; key = oligo coordinate,
        value = oligo sequence
    
    '''
    
    print('Loading reference fasta file...')
    seq_dict = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
    print('\t...complete\nGenerating oligos...')

    oligo_num = math.floor((max_dist-oligo)/step)
    oligo_seqs = {}
    assoc = {}
    
    with open(bed) as w:
        for x in w:
            chr_name, start, stop, name = x.strip().split('\t')
            
            chr_length = len(seq_dict[chr_name])
            start, stop = map(int, (start, stop))
            seq = seq_dict[chr_name].seq.upper()
            
            l_stop = start-10
            r_start = stop+10
        
            counter=1
            while counter<=oligo_num:
                
                l_start = l_stop-oligo
                l_coor = '{}:{}-{}'.format(chr_name, l_start, l_stop)
                r_stop = r_start+oligo
                r_coor = '{}:{}-{}'.format(chr_name, r_start, r_stop)
                
                if l_start<0:
                    print('Oligo {} for off-target site {} could not ' \
                          'be generated because it went beyond the start of ' \
                          'the chromosome'.format(l_coor, name),
                    file=sys.stderr)
                else:
                    l_seq = seq[l_start:l_stop]
                    oligo_seqs['{}-000-000-X'.format(l_coor)] = str(l_seq)
                    
                if r_stop>chr_length:
                    print('Oligo {} for off-target site {} could not ' \
                          'be generated because it went beyond the end of ' \
                          'the chromosome'.format(r_coor, name),
                    file=sys.stderr)
                else:
                    r_seq = seq[r_start:r_stop]
                    oligo_seqs['{}-000-000-X'.format(r_coor)] = str(r_seq)
                
                assoc[l_coor] = assoc[r_coor] = name
                
                l_stop-=step
                r_start+=step
                counter+=1
    
    pickle.dump(assoc, open('_tmp.p', 'wb'))
    print('\t...complete')
    return oligo_seqs

if __name__ == '__main__':
    
    import tools
    import argparse
    import subprocess
    
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
        '-b',
        '--bed',
        type = str,
        help = 'Path to bed file with off-target coordinates',
        required = True,
    )
    parser.add_argument(
        '-t',
        '--step_size',
        type = int,
        help = 'Step size (in bp) between adjacent oligos, default=10',
        default = 10,
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
        '-m',
        '--max_dist',
        type = int,
        help = 'The maximum distance away from the off-target site to ' \
               'design oligos to, default=200',
        default = 200,
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
    
    pass_seqs = gen_oligos(
        fa = args.fasta,
        bed = args.bed,
        step = args.step_size,
        oligo = args.oligo,
        max_dist = args.max_dist,
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
        subprocess.run('sort -k1,1 -k2,2n oligo_info.txt >all_oligos_info.txt',
                        shell=True)
        subprocess.run('rm -f oligo_info.txt', shell=True)
        subprocess.run('mv all_oligo_info.txt oligos_info.txt', shell=True)