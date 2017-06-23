#!/usr/bin/env python

import unittest
import os
import tiled
import regions
import re
import pybedtools
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import _verify_alphabet

#fa = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
#                  'Chromosomes/chr18.fa'
#full_chr = SeqIO.read(fa, 'fasta')
#self.start = 0
#self.stop = len(full_chr)

f = lambda x: re.split('\W+', x)
g = lambda x: tuple(map(int, re.split('\W+', x)[1:5]))

class OligoGenTest_Main(unittest.TestCase):
        
    def test_all_sequences_are_legitimate_dna_sequences(self):
        class DNAdict(): letters='GATCN'
        seq_list = [Seq(x, DNAdict) for x in self.seqs.values()]
        self.assertTrue(all(_verify_alphabet(seq) for seq in seq_list))
        
    def test_coordinates_can_be_coerced_to_int(self):
        full_list = [y for x in self.seqs for y in f(x)[1:5]]
        try:
            list(map(int, full_list))
        except ValueError:
            self.fail('Coordinates are not integers')
    
    def test_difference_between_oligo_start_stop_is_correct_length(self):
        self.assertTrue(
            all(g(x)[1]-g(x)[0]==self.oligo for x in self.seqs)
        )
        
    def test_all_oligos_are_correct_length(self):
        self.assertTrue(all(len(x) == self.oligo for x in self.seqs.values()))
        
    def test_sequence_is_same_if_using_bedtools_getfasta(self):
        a = pybedtools.BedTool(' '.join(f(next(iter(self.seqs)))[:3]), \
                               from_string=True)
        a = a.sequence(fi='/databank/igenomes/Mus_musculus/UCSC/mm9/' \
                          'Sequence/WholeGenomeFasta/genome.fa')
        seq_f = open(a.seqfn)
        self.assertMultiLineEqual(
            seq_f.read().split('\n')[1].upper(),
            next(iter(self.seqs.values()))
        )
        seq_f.close()
        
    def test_fasta_file_has_same_contents_as_oligo_dictionary(self):
        tiled.write_oligos(oligo_seqs=self.seqs)
        test_dict = {}
        with open('oligo_seqs.fa') as f:
            for x in f:
                test_dict[x.strip('>').rstrip('\n')] = next(f).rstrip('\n')
        os.remove('oligo_seqs.fa')
        self.assertDictEqual(self.seqs, test_dict)
        
class OligoGenTest_Capture(unittest.TestCase):
    
    def test_fragment_start_equals_left_oligo_start(self):
        self.assertTrue(
            all(g(x)[0]==g(x)[2] for x in self.seqs if f(x)[5]=='L'),
        )
    
    def test_fragment_stop_equals_right_oligo_stop(self):
        self.assertTrue(
            all(g(x)[1]==g(x)[3] for x in self.seqs if f(x)[5]=='R'),
        )
        
    def test_oligo_is_inside_fragment(self):
        self.assertTrue(
            all((g(x)[0]>=g(x)[2]) & (g(x)[1]<=g(x)[3]) for x in self.seqs),
        )
    
    def test_all_left_oligos_start_with_restriction_site(self):
        self.assertTrue(
            all(y[0:len(self.res_site)] == self.res_site \
                for x, y in self.seqs.items() if f(x)[5]=='L'),
        )
    
    def test_all_right_oligos_end_with_restriction_site(self):
        self.assertTrue(
            all(y[-len(self.res_site):] == self.res_site \
                for x, y in self.seqs.items() if f(x)[5]=='R'),
        )
        
class OligoGenTest_Blocks(unittest.TestCase):
    
    def test_all_oligo_start_coordinates_are_within_specified_range(self):
        self.assertTrue(
            all((g(x)[0] >= self.start) & (g(x)[0] <= self.stop) \
                for x in self.seqs),
        )
       
    def test_all_oligo_stop_coordinates_are_within_specified_range(self):
        self.assertTrue(
            all((g(x)[1] >= self.start) & (g(x)[1] <= self.stop) \
                for x in self.seqs),
        )
        
class OligoGenTest_Tiled(OligoGenTest_Main, OligoGenTest_Capture, OligoGenTest_Blocks):
    
    @classmethod   
    def setUpClass(self,
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
        self.seqs = tiled.gen_oligos_capture(
            fa = self.fa,
            chromosome = self.chromosome,
            enzyme = self.enzyme,
            oligo = self.oligo,
            region = self.region,
        )
    
    @classmethod
    def tearDownClass(self):
        pass

        
class OligoGenTest_FISH(OligoGenTest_Main, OligoGenTest_Blocks):
     
    @classmethod   
    def setUpClass(self,
              fa = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
                   'Chromosomes/chr18.fa',
              step=30,
              oligo=30,
              region = '44455000-44555000',
              chromosome = 18):
        self.fa = fa
        self.step = step
        self.oligo = oligo
        self.chromosome = chromosome
        self.region = region
        self.start, self.stop = map(int, region.split('-'))
        self.seqs = tiled.gen_oligos_fish(
            fa = self.fa,
            chromosome = self.chromosome,
            step = self.step,
            oligo = self.oligo,
            region = self.region,
        )
        
class OligoGenTest_Regions(OligoGenTest_Main, OligoGenTest_Capture):
    
    @classmethod   
    def setUpClass(self,
              fa = '/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
                   'WholeGenomeFasta/genome.fa',
              bed = '/t1-data1/WTSA_Dev/jkerry/CaptureC/DDownes/' \
                    'CapsequmInput_2.txt',
              enzyme = 'DpnII',
              oligo = 30):
        self.fa = fa
        self.bed = bed
        self.enzyme = enzyme
        self.oligo = oligo
        self.res_site = tiled.rs_dict[enzyme]
        self.seqs = regions.gen_oligos(
            fa = self.fa,
            bed = self.bed,
            enzyme = self.enzyme,
            oligo = self.oligo,
        )


if __name__ == '__main__':
    unittest.main()
    