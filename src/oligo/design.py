from __future__ import print_function, division

import argparse
import re
import sys

import numpy as np

from oligo.tools import Tools

recognition_seq = {'DpnII': 'GATC',
                   'NlaIII': 'CATG',
                   'HindIII': 'AAGCTT'}

def _check_value(values, labels):
    for value, label in zip(values, labels):
        if (value<1) | isinstance(value, float):
            raise ValueError(
                '{} must be an integer greater than 0'.format(label))
            break
    else:
        return None
    
def _compile_chr_regex(genome):
    if genome.startswith('hg'):
        twenties = '2[0-2]|'
    elif genome.startswith('mm'):
        twenties = ''
    valid_chroms = re.compile('^chr([1-9]|1[0-9]|'+twenties+'X|Y)$')
    
    return valid_chroms

def _create_key(chrom, *args):
    return (':'.join((chrom, '-'.join((map(str, x + ('000','000','X')))))) for x in args)

def _create_key_frag(chrom, oligo_coor, frag_coor, side):
    key = ':'.join((chrom, '-'.join(map(str, (oligo_coor + frag_coor)))))
    key = '-'.join((key, side))
            
    return key

def _get_sequence(seq, *args):
    return (str(seq[x[0]:x[1]]) for x in args)  # returns a generator so that all sequences outside of the dictionary are not retained in the memory

def _validate_chrom(chrom, regex):
    if not regex.match(chrom):
        raise ChromosomeError('Unrecognised chromosome {}. Skipping.')
    return None

class ChromosomeError(Exception):
    pass

class FragmentError(Exception):
    pass

class FragmentMixin(object):
    
    def _get_fragment_seqs(self, chrom, start, stop, chrom_seq=''):
        frag_length = stop - start
        if frag_length < self.oligo:
            raise FragmentError('{} is too small to design oligos in. Skipping.')
        
        left_coor = (start, start + self.oligo)
        right_coor = (stop - self.oligo, stop)
        frag_coor = (start, stop)
        
        left_key = _create_key_frag(chrom, left_coor, frag_coor, 'L')
        right_key = _create_key_frag(chrom, right_coor, frag_coor, 'R')
        
        if left_key in self.oligo_seqs:
            raise FragmentError('{} is redundant to another position. Skipping.')
        
        if not chrom_seq: chrom_seq = self.genome_seq[chrom].seq.upper()
        left_seq, right_seq = _get_sequence(chrom_seq, *(left_coor, right_coor))
        
        self.oligo_seqs[left_key] = left_seq     
        if frag_length > self.oligo: self.oligo_seqs[right_key] = right_seq
        
        return None

class Capture(Tools, FragmentMixin):
    """Designs oligos for Capture-C"""
    
    __doc__ += Tools.__doc__
        
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
            The enzyme for digestion, default = DpnII
        oligo : int, optional
            The length of the oligo to design (bp), default = 70
        
        Returns
        -------
        self : object
        
        """
        
        _check_value((oligo,), ('Oligo size',))
        self._create_attr(oligo)
        
        cut_sites = {}
        cut_size = len(recognition_seq[enzyme])
        rec_seq = re.compile(recognition_seq[enzyme])
        
        print('Generating oligos...')
        for chrom, attr in self.genome_seq.items():
            chrom_seq = str(attr.seq).upper()
            cut_sites[chrom] = [cut_site.start() for cut_site
                                in rec_seq.finditer(chrom_seq)]
            cut_sites[chrom] = np.array(cut_sites[chrom])
            
        valid_chroms = _compile_chr_regex(self.genome)
            
        with open(bed) as viewpoints:
            for vp in viewpoints:
                chrom, vp_start, vp_stop, name = vp.strip().split('\t')
                
                try:
                    _validate_chrom(chrom, valid_chroms)
                except ChromosomeError as e:
                    print(str(e).format(chrom), file=sys.stderr)
                    continue
                
                vp_start = int(vp_start)
                try:
                    frag_start = cut_sites[chrom][cut_sites[chrom] <= vp_start][-1]
                    frag_stop = cut_sites[chrom][cut_sites[chrom] > vp_start][0] + cut_size # currently this picks an adjacent fragment if the site is in a cutsite; are we okay with that?
                except IndexError:
                    print('Viewpoint {} is not in a closed fragment. '
                          'Skipping.'.format(name))
                    continue
                
                try:
                    self._get_fragment_seqs(chrom, frag_start, frag_stop)
                except FragmentError as e:
                    frag_id = '{}:{}-{} ({}) is in a fragment that'.format(
                        chrom, vp_start, vp_stop, name)
                    print(str(e).format(frag_id), file=sys.stderr)
                    continue
                
                frag_key = '{}:{}-{}'.format(chrom, frag_start, frag_stop)
                self._assoc[frag_key] = '{}{},'.format(
                    self._assoc.get(frag_key, ''), name)
        
        print('\t...complete.')
        if __name__ != '__main__':
            print('Oligos stored in the oligo_seqs attribute')
        
        return self
    
    def __str__(self):
        
        return 'Capture-C oligo design object for the {} genome'.format(
            self.genome)
    
class Tiled(Tools, FragmentMixin):
    """Designs oligos adjacent to each other or on adjacent fragments"""
    
    __doc__ += Tools.__doc__
    
    def gen_oligos_capture(self, chrom, region='', enzyme='DpnII', oligo=70):
        """Designs oligos for multiple adjacent restriction fragments
        across a specified region of a chromosome, or for the entire
        chromosome.
        
        Parameters
        ----------
        chrom : str
            Chromosome number/letter e.g. 7 or X
        region : str, optional
            The region of the chromosome to design oligos, e.g.
            10000-20000; omit this option to design oligos over the
            entire chromosome
        enzyme : {'DpnII', 'NlaIII', 'HindIII'}, optional
            The enzyme for digestion, default=DpnII
        oligo : int, optional
            The length of the oligos to design (bp), default=70
            
        Returns
        -------
        self : object
    
        """
        
        _check_value((oligo,), ('Oligo size',))
        self._create_attr(oligo)
        
        if 'chr' not in str(chrom): chrom = ''.join(('chr'+str(chrom)))
        chrom_seq = self.genome_seq[chrom].seq.upper()
            
        start, stop = (0, len(chrom_seq)) if not region else map(int, region.split('-'))
        
        rec_seq = re.compile(recognition_seq[enzyme])
        cut_size = len(recognition_seq[enzyme])
        cut_sites = [cut_site.start()+start for cut_site
                     in rec_seq.finditer(str(chrom_seq[start:stop]))]
        
        print('Generating oligos...')   
        for i in range(len(cut_sites)-1):
            j = i + 1
            
            frag_start = cut_sites[i]
            frag_stop = cut_sites[j]+cut_size 
            
            try:
                self._get_fragment_seqs(chrom, frag_start, frag_stop, chrom_seq=chrom_seq)
            except FragmentError as e:
                frag_id = 'The fragment {}:{}-{}'.format(chrom, frag_start, frag_stop)
                print(str(e).format(frag_id), file=sys.stderr)
                continue
                
        print('\t...complete.')
        if __name__ != '__main__':
            print('Oligos stored in the oligo_seqs attribute')
        
        return self
    
    def gen_oligos_contig(self, chrom, region='', step=70, oligo=70):
        """Designs adjacent oligos based on a user-defined step size,
        across a specified region of a chromosome, or for the entire
        chromosome.
        
        Parameters
        ----------
        chrom : str
            Chromosome number/letter e.g. 7 or X
        region : str, optional
            The region of the chromosome to design oligos, e.g.
            10000-20000; omit this option to design oligos over the
        step : int, optional
            The step size, or how close adjacent oligos should be, default=70;
            for end-to-end oligos, set `step` to equal `oligo`
        oligo : int, optional
            The length of the oligos to design (bp), default=70
            
        Raises
        ------
        ValueError
            If `step` or `oligo` <1 or not an integer
            
        Returns
        -------
        self : object
        
        """
        
        _check_value((step, oligo), ('Step size', 'Oligo size'))
        self._create_attr(oligo)
        
        if not chrom.startswith('chr'): chrom = ''.join(('chr'+str(chrom)))
        chrom_seq = self.genome_seq[chrom].seq.upper()
            
        start, stop = (0, len(chrom_seq)) if not region else map(
            int, region.split('-'))
        stop = stop - oligo
        
        print('Generating oligos...')
        coors = [(x, x + oligo) for x in range(start, stop+1, step)]
        sequences = _get_sequence(chrom_seq, *coors)
        keys = _create_key(chrom, *coors)
        self.oligo_seqs.update(zip(keys, sequences))
        
        print('\t...complete.')
        if __name__ != '__main__':
            print('Oligos stored in the oligo_seqs attribute')
        
        return self
    
    def __str__(self):
        
        return 'Tiled Capture oligo design object for the {} genome'.format(
            self.genome)
    
class OffTarget(Tools):
    """Designs oligos adjacent to potential CRISPR off-target sites"""
    
    __doc__ += Tools.__doc__
    
    def gen_oligos(self, bed, oligo=70, step=10, max_dist=200):
        r"""Designs oligos adjacent to user-supplied coordinates for
        potential CRISPR off-target cleavage sites
        
        Parameters
        ----------
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
            oligos to, default = 200
        
        Returns
        -------
        self : object
        
        """
        
        _check_value((step, oligo, max_dist),
                    ('Step size', 'Oligo size', 'Maximum distance'))
        self._create_attr(oligo)
    
        valid_chroms = _compile_chr_regex(self.genome)
        print('Generating oligos...')
        with open(bed) as crispr_sites:
            for site in crispr_sites:
                chrom, start, stop, name = site.strip().split('\t')
                
                try:
                    _validate_chrom(chrom, valid_chroms)
                except ChromosomeError as e:
                    print(str(e).format(chrom), file=sys.stderr)
                    continue
                  
                start, stop = map(int, (start, stop))
                chrom_seq = self.genome_seq[chrom].seq.upper()
                chrom_length = len(chrom_seq)
                
                exc_start = start-10
                exc_stop = stop+10
                
                upstream = [(x-oligo, x) for x in
                    range(start-max_dist+oligo, exc_start+1, step)
                    if x-oligo>=0]
                downstream = [(x, x+oligo) for x in
                    range(exc_stop, stop+max_dist-oligo+1, step)
                    if x+oligo<=chrom_length]
                all_coors = upstream+downstream
                sequences = _get_sequence(chrom_seq, *all_coors)
                keys = list(_create_key(chrom, *all_coors))
                self.oligo_seqs.update(zip(keys, sequences))
                self._assoc.update({x[:-10]: name for x in keys})
        
        print('\t...complete.')
        if __name__ != '__main__':
            print('Oligos stored in the oligo_seqs attribute')
        
        return self
    
    def __str__(self):
        
        return 'OffTarget Capture oligo design object for the {} genome'.format(
            self.genome)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    args = parser.parse_args(sys.argv[2:])
    
    if not args.blat and not args.star_index:
        msg = '-s/--star_index argument is required if --blat is not selected'
        parser.error(msg)
