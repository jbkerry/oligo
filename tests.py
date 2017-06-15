#!/usr/bin/env python

import unittest
import os
import tiled
import re
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import _verify_alphabet

#fa = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
#                  'Chromosomes/chr18.fa'
#full_chr = SeqIO.read(fa, 'fasta')
#self.start = 0
#self.stop = len(full_chr)

class OligoGenTest(unittest.TestCase):
    
    def setUp(self,
              fa = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
                   'Chromosomes/chr18.fa',
              enzyme='DpnII',
              oligo=30,
              region = '44455000-44555000',
              chromosome = 18):
        self.fa = fa
        self.enzyme = enzyme
        self.res_site = tiled.rs_dict[enzyme]
        self.oligo = oligo
        self.chromosome = chromosome
        self.region = region
        self.start, self.stop = map(int, region.split('-'))
        
        tiled.gen_oligos_capture(
            fa = self.fa,
            chromosome = self.chromosome,
            enzyme = self.enzyme,
            oligo = self.oligo,
            region = self.region,
        )
        with open('oligo_seqs.fa') as f:
            self.lines = [x.rstrip('\n') for x in f]
    
    def tearDown(self):
        os.remove('oligo_seqs.fa')
    
    def test_oligo_fasta_created(self):
        self.assertTrue(os.path.exists('oligo_seqs.fa'))
        
    def test_all_odd_numbered_lines_are_fasta_name_lines(self):
        name_list = self.lines[::2]
        self.assertTrue(all(x[0] == '>' for x in name_list))
        
    def test_all_even_numbered_lines_are_dna_sequences(self):
        class DNAdict(): letters='GATCN'
        seq_list = [Seq(x, DNAdict) for x in self.lines[1::2]]
        self.assertTrue(all(_verify_alphabet(seq) for seq in seq_list))
        
    def test_fasta_file_has_an_even_number_of_lines(self):
        self.assertEqual(len(self.lines) % 2, 0)
        
    def test_coordinates_can_be_coerced_to_int(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        full_list = [y for x in coor_list for y in re.split('\W+', x)[1:5]]
        #full_list = []
        #for x in coor_list: full_list += re.split('\W+', x)[1:5]
        try:
            list(map(int, full_list))
        except ValueError:
            self.fail('Coordinates are not integers')
    
    def test_difference_between_oligo_start_stop_is_correct_length(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        oligo_starts = [re.split('\W+', x)[1] for x in coor_list]
        oligo_stops = [re.split('\W+', x)[2] for x in coor_list]
        self.assertTrue(
            all(int(y)-int(x)==self.oligo \
               for x, y in zip(oligo_starts, oligo_stops))
        )
    
    def test_fragment_start_equals_oligo_start(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        oligo_starts = [re.split('\W+', x)[1] for x in coor_list \
                        if re.split('\W+', x)[5]=='L']
        frag_starts = [re.split('\W+', x)[3] for x in coor_list \
                       if re.split('\W+', x)[5]=='L']
        self.assertListEqual(oligo_starts, frag_starts)
    
    def test_fragment_stop_equals_oligo_stop(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        oligo_stops = [re.split('\W+', x)[2] for x in coor_list \
                        if re.split('\W+', x)[5]=='R']
        frag_stops = [re.split('\W+', x)[4] for x in coor_list \
                        if re.split('\W+', x)[5]=='R']
        self.assertListEqual(oligo_stops, frag_stops)
        
    def test_oligo_is_inside_fragment(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        full_list = [list(map(int, re.split('\W+', x)[1:5])) \
                     for x in coor_list]
        self.assertTrue(
            all((x[0]>=x[2]) & (x[1]<=x[3]) for x in full_list)
        )
        
    def test_there_are_no_duplicate_coordinates(self):
        name_list = self.lines[::2]
        self.assertEqual(len(name_list), len(set(name_list)))
        
    def test_all_oligos_are_correct_length(self):
        seq_list = self.lines[1::2]
        self.assertTrue(all(len(x) == self.oligo for x in seq_list))
        
    def test_all_oligo_start_coordinates_are_within_specified_range(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        oligo_start_list = [int(re.split('\W+', x)[1]) for x in coor_list]
        self.assertTrue(
            all((x >= self.start) & (x <= self.stop) for x in oligo_start_list),
        )
        
    def test_all_oligo_stop_coordinates_are_within_specified_range(self):
        coor_list = [x.strip('>') for x in self.lines[::2]]
        oligo_stop_list = [int(re.split('\W+', x)[2]) for x in coor_list]
        self.assertTrue(
            all((x >= self.start) & (x <= self.stop) for x in oligo_stop_list),
        )
        
    def test_all_left_oligos_start_with_restriction_site(self):
        left_seqs = []; i = 0
        for x in self.lines[::2]:
            if re.split('\W+', x)[-1] == 'L': left_seqs.append(self.lines[i+1])
            i+=2
        self.assertTrue(
            all(x[0:len(self.res_site)] == self.res_site for x in left_seqs),
        )
    
    def test_all_right_oligos_end_with_restriction_site(self):
        right_seqs = []; i = 0
        for x in self.lines[::2]:
            if re.split('\W+', x)[-1] == 'R': right_seqs.append(self.lines[i+1])
            i+=2
        self.assertTrue(
            all(x[-len(self.res_site):] == self.res_site for x in right_seqs),
        )
        
    def test_sequence_is_same_if_using_bedtools_getfasta(self):
        first_coor, first_seq = self.lines[:2]
        coor_list = re.split('\W+', first_coor.strip('>'))
        a = pybedtools.BedTool(' '.join(coor_list[:3]), from_string=True)
        a = a.sequence(fi='/databank/igenomes/Mus_musculus/UCSC/mm9/' \
                          'Sequence/WholeGenomeFasta/genome.fa')
        f = open(a.seqfn)
        self.assertMultiLineEqual(
            f.read().split('\n')[1].upper(),
            first_seq
        )
        f.close()

if __name__ == '__main__':
    unittest.main()
    