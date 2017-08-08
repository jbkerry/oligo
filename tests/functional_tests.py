#!/usr/bin/env python

import os
import unittest

import tiled
from Bio.Seq import Seq
from Bio.Alphabet import _verify_alphabet

class OligoGenTest(unittest.TestCase):
    
    @classmethod   
    def setUpClass(self):
        os.system('python tiled.py -f /databank/igenomes/Mus_musculus/UCSC/' \
                  'mm9/Sequence/Chromosomes/chr10.fa -g mm9 -c 10 ' \
                  '-r 12400000-12450000 -s /databank/igenomes/Mus_musculus/UCSC/' \
                  'mm9/Sequence/STAR --test_fasta')
        
        with open('oligo_seqs.fa') as f:
            self.lines = [x.rstrip('\n') for x in f]
    
    @classmethod
    def tearDownClass(self):
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
        
    def test_there_are_no_duplicate_coordinates(self):
        name_list = self.lines[::2]
        self.assertEqual(len(name_list), len(set(name_list)))
        
if __name__ == '__main__':
    unittest.main()
    