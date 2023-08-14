from __future__ import print_function, division

from collections import namedtuple
import os
import re
import subprocess

import pandas as pd
import pysam
from Bio import SeqIO


GENOME_MAP = {
    "mm": {"species": "Mus musculus", "docker_lib_file": "/usr/local/mouse.hmm"},
    "hg": {"species": "Homo sapiens", "docker_lib_file": "/usr/local/humans.hmm"},
}

blat_param = '-stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999'
star_param = '--readFilesIn {} --genomeDir {} --runThreadN 4 --genomeLoad NoSharedMemory ' \
             '--outFilterMultimapScoreRange 1000 --outFilterMultimapNmax ' \
             '100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 ' \
             '--seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax ' \
             '15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 ' \
             '--limitOutSAMoneReadBytes 400000 --outFileNamePrefix oligos_'

pat = re.compile('^[A-Z]')

class Tools(object):
    """
    
    Parameters
    ----------
    genome : {'mm9', 'mm10', 'hg18', 'hg19', 'hg38'}
        Genome build
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
        Contains all oligo sequences after generating oligos
        
    """
    
    def __init__(self, genome, fa, config_path, blat=False):
        self.genome = genome
        self.fa = fa
        self.paths = dict((x.strip().split(' = ') for x in open(config_path) if pat.match(x)))
        self.blat = blat
        self.fasta = 'oligo_seqs.fa'
        if self.__class__.__name__ != 'Tools':
            print('Loading reference fasta file...')
            self.genome_seq = SeqIO.to_dict(SeqIO.parse(fa, 'fasta'))
            print('\t...complete')
            
    def _create_attr(self, oligo):
        """Creates `oligo`, `oligo_seqs` and `_assoc` attributes"""
        
        self.oligo = oligo
        self.oligo_seqs = {}
        self._assoc = {}

    def write_fasta(self):
        """Writes `oligo_seqs` attribute to fasta file"""
        
        with open(self.fasta, 'w') as fa_w:
            for key, value in self.oligo_seqs.items():
                fa_w.write('>{}\n{}\n'.format(key, value))
        
        print('Wrote oligos to {}'.format(self.fasta))
        
        return None
    
    def detect_repeats(self):
        """Detects repeat sequences in oligos, using RepeatMasker"""
        
        options = ('RM_PATH', 'RepeatMasker', 'RepeatMasker',
                   'rm_log.txt', ''.join((self.fasta, '.out')))
        genome_id = self.genome.lower()[:2]
        cmd = f"-noint -s -species {GENOME_MAP[genome_id]['species']} "
        if os.getenv("USE_CUSTOM_RM_LIB"):
            cmd += f"-lib {GENOME_MAP[genome_id]['docker_lib_file']} "
        cmd += f"{self.fasta}"
        msg = 'Checking for repeat sequences in oligos,'
        
        self._run_command(options, cmd, msg)
        
        return self
    
    def align_to_genome(self, s_idx=''):
        """Aligns oligos to the genome using BLAT or STAR
        
        Parameters
        ----------
        s_idx : str
            Path to the directory containing the STAR index for this
            genome (not required if blat=True)
            
        Raises
        ------
        AttributeError
            If `blat` = False but `s_idx` is not specified
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
        
        if self.blat:
            blat_out = 'blat_out.psl'
            options = ('BLAT_PATH', 'blat', 'BLAT', 'blat_log.txt', blat_out)
            cmd = ' '.join((blat_param, self.fa, self.fasta, blat_out))
        else:
            options = ('STAR_PATH', 'STAR', 'STAR', 'star_log.txt',
                       'oligos_Aligned.out.sam')
            cmd = star_param.format(self.fasta, s_idx)
        msg = 'Aligning oligos to the genome,'
        
        self._run_command(options, cmd, msg)
        
        return None
    
    def extract_repeats(self):
        """Extracts information of repeat content from RepeatMasker output
        file for every oligo
        
        """
        
        try:
            self._oligo_stats
        except AttributeError:
            self._populate_oligo_stats()
        
        with open('.'.join((self.fasta, 'out'))) as repeats_file:
            if len(repeats_file.readlines())>1:
                repeats_file.seek(0)
                for _ in range(3):
                    next(repeats_file)
                for line in repeats_file:
                    parts = re.split("\s+", line.strip())
                    oligo_name = parts[4]
                    repeat_type = parts[9]
                    fragment_side = re.split("\W+", oligo_name)[5]
                    if len(fragment_side)>1:
                        oligo_name = oligo_name.split('_')[0]
                    
                    qstart, qstop = map(int, (parts[5:7]))
                    length = (qstop-qstart) + 1
                    if length > self._oligo_stats[oligo_name]['repeat_length']:
                        self._oligo_stats[oligo_name]['repeat_length'] = length
                        self._oligo_stats[oligo_name]['repeat_type'] = repeat_type
                msg = 'Repeat scores calculated'
            else:
                msg = 'No repeats detected'
        
        print(msg)
        
        return self
    
    def calculate_density(self,
                          sam='oligos_Aligned.out.sam',
                          blat_file='blat_out.psl'):
        """Calculates the repeat scores and off-target binding for
        each oligo based on their scores from RepeatMasker and
        STAR/BLAT. Outputs results to `oligo_info.txt`.
        
        Parameters
        ----------
        sam : str
            Path to STAR alignment (.sam) file from `align_to_genome`
            (not required if `blat`=True), default = oligos_Aligned.out.sam
        blat_file : str
            Path to BLAT alignment (.psl) file from `align_to_genome`
            (not required if `blat`=False), default = blat_out.psl
        
        """
        
        try:
            self._oligo_stats
        except AttributeError:
            self._populate_oligo_stats()
        
        if self.blat:        
            with open(blat_file) as f:
                for _ in range(5):
                    next(f)
                for line in f:
                    parts = re.split("\s+", line.strip())
                    oligo_name = parts[9]
                    qgapbases, qstart, qend = map(int, (parts[5], parts[11], 
                                                        parts[12]))
                    self._oligo_stats[oligo_name]['multimap'] += 1
                    self._oligo_stats[oligo_name]['matches'] += (int(qend) -
                                                           int(qstart)) + 1
                    self._oligo_stats[oligo_name]['mismatches'] += int(qgapbases)  
        else:
            sf = pysam.AlignmentFile(sam, 'r')
            for r in sf.fetch(until_eof=True):
                oligo_name = r.query_name
                if self._oligo_stats[oligo_name]['multimap'] == 0:
                    self._oligo_stats[oligo_name]['multimap'] = r.get_tag('NH')
                        
                for block in r.cigartuples:
                    if block[0] == 0:
                        self._oligo_stats[oligo_name]['matches'] += block[1]
                    elif (block[0] == 1) | (block[0] == 2):
                        self._oligo_stats[oligo_name]['mismatches'] += block[1]
            
        for oligo in self._oligo_stats:
            score = self._oligo_stats[oligo]['matches'] - self._oligo_stats[oligo]['mismatches']
            density = score / len(self._oligo_stats[oligo]['sequence'])
            self._oligo_stats[oligo]['density'] = float("{0:.2f}".format(density))      
        
        print('Density scores calculated')
        
        return self
    
    def write_oligo_info(self):
        """Writes oligo stats to oligo_info.txt"""
        
        p = re.compile('\W+')
        with open('oligo_info.txt', 'w') as output:
            output.write('chr\tstart\tstop\tfragment_start\tfragment_stop\t'
                    'side_of_fragment\tsequence\ttotal_number_of_alignments\t'
                    'density_score\trepeat_length\trepeat_class\tGC%\t'
                    'associations\n')
            for oligo, stats in self._oligo_stats.items():
                oligo_parts = (chrom, read_start, read_stop, frag_start,
                               frag_stop, frag_side) = p.split(oligo)
                
                has_fragment = True
                if (frag_start, frag_stop, frag_side) == ('000', '000', 'X'):
                    has_fragment = False
                    #frag_start, frag_stop, frag_side = '.' * 3
                    oligo_parts[3:] = '.' * 3
                
                try:
                    self._assoc
                except AttributeError:
                    associations = '.'
                else:
                    if self._assoc:
                        if has_fragment:
                            coor = '{}:{}-{}'.format(chrom, frag_start, frag_stop)
                        else:
                            coor = '{}:{}-{}'.format(chrom, read_start, read_stop)
                        associations = self._assoc.get(coor, '.')
                    else:
                        associations = '.'
                
                keys = ('sequence', 'multimap', 'density', 'repeat_length',
                        'repeat_type', 'GC%')
                to_write = oligo_parts + [str(stats[x])
                                          for x in keys] + [associations]
                output.write('{}\n'.format('\t'.join(to_write)))
    
        sorted_df = self._sort_file()
        sorted_df.to_csv('oligo_info.txt', sep='\t', index=False, na_rep='NA')
        print('Oligo information written to oligo_info.txt')
        
        return None
    
    def _run_command(self, options, cmd, msg):
        """Runs a command using subprocess"""
        
        CmdOptions = namedtuple('CmdOptions', ['paths_key', 'exe', 'name',
                                               'log_file', 'output_file'])           
        run_options = CmdOptions._make(options)
        path = os.path.join(self.paths[run_options.paths_key], run_options.exe)
        print('{} with {}...'.format(msg, run_options.name))
        log = open(run_options.log_file, 'w')
        subprocess.call(' '.join((path, cmd)), shell=True, stdout=log,
                        stderr=log)
        log.close()
        print('\t...complete. Output written to {}'.format(
            run_options.output_file))
        
        return None
    
    def _populate_oligo_stats(self):
        """Populates _oligo_stats attribute with default values""" 
        self._oligo_stats = {}
        with open(self.fasta) as fasta_file:
            for line in fasta_file:
                oligo_name = line.lstrip('>').strip()
                read_seq = next(fasta_file).strip()
                self._oligo_stats[oligo_name] = {
                    'sequence': read_seq,
                    'multimap': 0,
                    'density': 0,
                    'repeat_length': 0,
                    'repeat_type': 'NA',
                    'GC%': self._get_gc(read_seq),
                    'matches': 0,
                    'mismatches': 0,
                }
        
        return None
        
    def _get_gc(self, x):
        """Calculates GC percentage of a DNA sequence"""
        
        gc_decimal = (x.count('C') + x.count('G'))/len(x)
        gc_decimal = float("{0:.2f}".format(gc_decimal))
        gc_perc = int(gc_decimal*100)
        
        return gc_perc
    
    def _sort_file(self):
        """Sorts oligo output file"""
        df = pd.read_table('oligo_info.txt', header=0)
        #df['chr'] = [int(x[3:]) for x in df['chr']]  # this threw an error for chromosomes X and Y
        df.sort_values(['chr', 'start'], inplace=True)
        #df['chr'] = 'chr' + df['chr'].map(str)
        
        return df
    
    def __repr__(self):
        
        return '{}(genome={}, fa={}, blat={})'.format(self.__class__.__name__,
                                                      self.genome,
                                                      self.fa,
                                                      self.blat)
