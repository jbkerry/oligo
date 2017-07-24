#!/usr/bin/env python

import os
import subprocess
import re
import sys
import pickle

import pysam
import numpy as np
from Bio import SeqIO

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

spe_dict = {'mm9': 'mouse',
            'mm10': 'mouse',
            'hg18': 'human',
            'hg19': 'human',
            'hg38': 'human'}

blat_param = '-stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999'
star_param = '--runThreadN 4 --genomeLoad NoSharedMemory ' \
             '--outFilterMultimapScoreRange 1000 --outFilterMultimapNmax ' \
             '100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 ' \
             '--seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax ' \
             '15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 ' \
             '--limitOutSAMoneReadBytes 400000 --outFileNamePrefix tiled_'

class UtilsMixin(object):
    """Contains a suite of methods inherited by Capture, Tiled and OffTarget"""

    def write_oligos(self):
        """Writes `oligo_seqs` attribute to fasta file"""
        
        with open(self.fasta, 'w') as fa_w:
            for key, value in self.oligo_seqs.items():
                fa_w.write('>{}\n{}\n'.format(key, value))
        
        print('Wrote oligos to {}'.format(self.fasta))
        
        return None
    
    def check_off_target(self, s_idx=''):
        """Checks for repeat sequences in oligos in fasta file using
        RepeatMasker and checks for off-target binding using either
        BLAT or STAR
        
        Parameters
        ----------
        s_idx : str
            Path to the directory containing the STAR index for this
            genome (not required if blat=True)
            
        Raises
        ------
        AttributeError
            If `blat`=False but `s_idx` is not specified
        FileNotFoundError
            If a fasta file with the name specified by the `fasta`
            attribute is not found
            
        """
        
        if (not self.blat) and (not s_idx):
            raise AttributeError('Path to STAR index must be set if '
                                 'blat=False')
        if not os.path.exists(self.fasta):
            raise FileNotFoundError('A valid FASTA file with the name {} was '
                                    'not found'.format(self.fasta))
        
        p = re.compile('^[A-Z]')
        config_file = './config.txt'
        path_list = [x.rstrip('\n') for x in open(config_file) if p.match(x)]
        path_dict = dict(item.split(' = ') for item in path_list)
        rm_path = os.path.join(path_dict['RM_PATH'], 'RepeatMasker')
        print('Checking for repeat sequences in oligos...')
        rm_out = open('rm_log.txt', 'w')
        subprocess.run(
            '{} -noint -s -species {} {}'.format(rm_path,
                                                 spe_dict[self.genome.lower()],
                                                 self.fasta),
            shell = True,
            stdout = rm_out,
            stderr = rm_out,
        )
        rm_out.close()
        print('\t...complete')
        
        if self.blat:
            path = os.path.join(path_dict['BLAT_PATH'], 'blat')
            print('Checking off-target binding with BLAT...')
            blat_out = open('blat_log.txt', 'w')
            subprocess.run(
                '{} {} {} {} blat_out.psl'.format(path, blat_param,
                                                  self.fa, self.fasta),
                shell=True,
                stdout = blat_out,
                stderr = blat_out,
            )
            blat_out.close()
        else:
            path = os.path.join(path_dict['STAR_PATH'], 'STAR')
            print('Checking off-target binding with STAR...')
            star_out = open('star_log.txt', 'w')
            subprocess.run(
                '{} --readFilesIn {} --genomeDir {} {}'.format(path,
                                                               self.fasta,
                                                               s_idx,
                                                               star_param),
                shell = True,
                stdout = star_out,
                stderr = star_out,
            )
            star_out.close()
            
        print('\t...complete')
        
        return None
    
    def get_density(self,
                    sam='tiled_Aligned.out.sam',
                    blat_file='blat_out.psl'):
        """Calculates the and repeat scores and off-target binding for
        each oligo based on their scores from RepeatMasker and
        STAR/BLAT. Outputs results to `oligo_info.txt`.
        
        Parameters
        ----------
        sam : str
            Path to STAR alignment (.sam) file from `check_off_targets`
            (not required if `blat`=True), default = tiled_Aligned.out.sam
        blat_file : str
            Path to BLAT alignment (.psl) file from `check_off_targets`
            (not required if `blat`=False), default = blat_out.psl
        
        """
        
        #all_oligos[oligo] = [sequence, nh, density, rep_len,
        #                     rep_type, gc_perc, matches, mismatches]
        self._all_oligos = {}
        with open(self.fasta) as f:
            for x in f:
                header = re.sub('>', '', x.rstrip('\n'))
                seq = next(f).rstrip('\n')
                self._all_oligos[header] = [seq, 0, 0, 0, 'NA', self._get_gc(x=seq), 0, 0]
        if self.blat:        
            with open(blat_file) as f:
                for _ in range(5):
                    next(f)
                for x in f:
                    parts = re.split("\s+", x.rstrip('\n'))
                    query = parts[9]
                    qgapbases, qstart, qend = map(int, (parts[5], parts[11], 
                                                           parts[12]))
                    self._all_oligos[query][1]+=1
                    self._all_oligos[query][6]+=(int(qend)-int(qstart))+1
                    self._all_oligos[query][7]+=int(qgapbases)  
        else:
            sf = pysam.AlignmentFile(sam, 'r')
            for r in sf.fetch(until_eof=True):
                if self._all_oligos[r.query_name][1] == 0:
                    self._all_oligos[r.query_name][1] = r.get_tag('NH')
                        
                for block in r.cigartuples:
                    if block[0]==0:
                        self._all_oligos[r.query_name][6]+=block[1]
                    elif (block[0]==1) | (block[0]==2):
                        self._all_oligos[r.query_name][7]+=block[1]
            
        for o in self._all_oligos:
            score = self._all_oligos[o][6]-self._all_oligos[o][7]
            density = score/len(self._all_oligos[o][0])
            self._all_oligos[o][2] = float("{0:.2f}".format(density))      
        
        print('Density scores calculated')
        self._get_repeats()
        self._write_file()
        
        return None
        
    def _get_gc(self, x):
        """Calculates GC percentage of a DNA sequence"""
        
        gc_perc = (x.count('C') + x.count('G'))/len(x)
        gc_perc = float("{0:.2f}".format(gc_perc))
        
        return gc_perc
        
    def _get_repeats(self):
        """Extracts information of repeat content from RepeatMasker output
        file for every oligo
        
        """
        
        with open(self.fasta+'.out') as f:
            if len(f.readlines())>1:
                f.seek(0)
                for _ in range(3):
                    next(f)
                for x in f:
                    parts = re.split("\s+", x.rstrip('\n'))
                    qname = parts[5]
                    rep_type = parts[10]
                    chr_name, start, stop, fragstart, fragend, side = re.split(
                                                                "\W+", qname)
                    if len(side)>1:
                        qname, dup = qname.split('_')
                    
                    qstart, qstop = map(int, (parts[6:8]))
                    length = (qstop - qstart)+1
                    if length>self._all_oligos[qname][3]:
                        self._all_oligos[qname][3:5] = length, rep_type
                msg = 'Repeat scores calculated'
            else:
                msg = 'No repeats detected'
        
        print(msg)
        
        return None
    
    def _write_file(self):
        """Sorts and writes oligo information to oligo_info.txt"""
        
        p = re.compile('\W+')
        assoc = ''
        if os.path.exists('_tmp.p'): assoc = pickle.load(open('_tmp.p', 'rb'))
        with open('oligo_info.txt', 'w') as f:
            f.write('chr\tstart\tstop\tfragment_start\tfragment_stop\t'
                    'side_of_fragment\tsequence\ttotal_number_of_alignments\t'
                    'density_score\trepeat_length\trepeat_class\tGC%\t'
                    'associations\n')
            for key, idx in self._all_oligos.items():
                w_list = p.split(key)+idx[:-2]
                if w_list[3:6]==['000', '000', 'X']: w_list[3:6] = '.'*3
                if assoc:
                    if w_list[3:6]==['000', '000', 'X']:
                        coor = '{}:{}-{}'.format(w_list[0], w_list[1], w_list[2])
                    else:
                        coor = '{}:{}-{}'.format(w_list[0], w_list[3], w_list[4])
                    w_list.append(assoc.get(coor, '.'))
                else:
                    w_list.append('.')
                f.write('\t'.join(map(str, w_list))+'\n')
        
        print('Oligo information written to oligo_info.txt')
        
        return None

class Capture(UtilsMixin):
    """Generates oligos for Capture-C
    
    Parameters
    ----------
    genome : str
        Genome build e.g. mm10 or hg38
    fa : str
        Path to reference genome fasta
    
    Attributes
    ----------
    blat : bool
        Check off-target binding using BLAT instead of STAR (not
        recommended for large designs), default = False
    fasta : str
        Name of fasta file for oligo sequences, default = oligo_seqs.fa
    oligo_seqs : dict
        Contains all oligo sequences after generating oligos; key =
        oligo coordinate, value = oligo sequence
        
    """
    
    def __init__(self, genome, fa, blat=False):
        self.genome = genome
        self.fa = fa
        self.blat = blat
        self.fasta = 'oligo_seqs.fa'
        
    def gen_oligos(self, bed, enzyme='DpnII', oligo=70):
        r"""Generates oligos flanking restriction fragments that encompass the
        coordinates supplied in the bed file
        
        Parameters
        ----------
        bed : str
            Path to tab-delimited bed file containing a list of coordinates for
            viewpoints in the capture experiment. Must be in the format
            'chr'\\t'start'\\t'stop'\\t'name'\\n
        enzyme : {'DpnII', 'NlaIII', 'HindIII'}, optional
            The enzyme for digestion, default=DpnII
        oligo : int, optional
            The length of the oligo to design (bp), default=70
        
        Returns
        -------
        self : object
        
        """
        
        cut_sites = {}
        self.oligo_seqs = {}
        assoc = {}
        cut_size = len(rs_dict[enzyme])
        p = re.compile(rs_dict[enzyme])
        
        print('Loading reference fasta file...')
        seq_dict = SeqIO.to_dict(SeqIO.parse(self.fa, 'fasta'))
        for x in seq_dict:
            seq = str(seq_dict[x].seq).upper()
            cut_sites[x] = [m.start() for m in p.finditer(seq)]
            cut_sites[x] = np.array(cut_sites[x])
        print('\t...complete\nGenerating oligos...')
        
        with open(bed) as w:
            for x in w:
                chr_name, start, stop, name = x.strip().split('\t')
                
                if '_' in chr_name: continue
                
                vp_coor = '{}:{}-{}'.format(chr_name, start, stop)
                start, stop = map(int, (start, stop))
                seq = seq_dict[chr_name].seq.upper()
                
                l_start = cut_sites[chr_name][cut_sites[chr_name]<=start][-1]
                r_stop = cut_sites[chr_name][cut_sites[chr_name]>=start][0] + cut_size # currently this picks an adjacent fragment if the site is in a cutsite; are we okay with that?
                frag_len = r_stop - l_start
                
                if frag_len>=oligo:
                    l_stop = l_start + oligo
                    r_start = r_stop - oligo
                    
                    l_tup = (l_start, l_stop, l_start, r_stop)
                    r_tup = (r_start, r_stop, l_start, r_stop)
                    l_seq = seq[l_start:l_stop]
                    r_seq = seq[r_start:r_stop]
                    
                    frag_key = '{}:{}-{}'.format(chr_name, l_start, r_stop)
                    assoc[frag_key] = '{}{},'.format(assoc.get(frag_key, ''), name)
                    
                    l_key = '{}:{}-L'.format(chr_name, '-'.join(map(str, l_tup)))
                    if l_key in self.oligo_seqs:
                        print('{} ({}) is redundant to another position'.format(
                            vp_coor, name), file=sys.stderr)
                        #not_done[vp_coor] = (name, 0) # 0 = redundant
                        continue
                    
                    self.oligo_seqs[l_key] = str(l_seq)
                    if frag_len>oligo:
                        self.oligo_seqs['{}:{}-R'.format(chr_name,
                                               '-'.join(map(str, r_tup)))
                                ] = str(r_seq)
                else:
                    print('{} ({}) was in a fragment that was too small'.format(
                            vp_coor, name), file=sys.stderr)
        
        pickle.dump(assoc, open('_tmp.p', 'wb'))
        print('\t...complete')
        
        return self
    
if __name__ == '__main__':
    try:
        class_arg = sys.argv[1]
    except IndexError:
        raise('oligo.py must be followed by a class name: Capture, Tiled or OffTarget')
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
        '-o',
        '--oligo',
        type = int,
        help = 'The size (in bp) of the oligo to design, default=70',
        default = 70,
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
    
    if class_arg == 'Capture': 
        parser.add_argument(
            '-b',
            '--bed',
            type = str,
            help = 'Path to bed file with capture viewpoint coordinates',
            required = True,
        )
        parser.add_argument(
            '-e',
            '--enzyme',
            type = str,
            help = 'Name of restriction enzyme, default=DpnII',
            default = 'DpnII',
            required = False,
        )
        
        args = parser.parse_args()
        c = Capture(genome=args.genome, fa=args.fasta)
        c = c.gen_oligos(
            bed = args.bed,
            enzyme = args.enzyme,
            oligo = args.oligo,
        )
    
    c.write_oligos()
    c.check_off_target(s_idx=args.star_index)
    c.get_density()
    

