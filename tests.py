#!/usr/bin/env python

import unittest
import os
import tiled
import re

class OligoGenTest(unittest.TestCase):
    
    def setUp(self):
        self.oligo = 30
        self.start = 44455000
        self.stop = 44555000
        tiled.gen_oligos_capture(
            fa='/databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/' \
               'Chromosomes/chr18.fa',
            chromosome = 18,
            oligo = self.oligo,
            region = '{}-{}'.format(self.start, self.stop)
        )
        with open('oligo_seqs.fa') as f:
            self.lines = [x.rstrip('\n') for x in f]
    
    def tearDown(self):
        os.remove('oligo_seqs.fa')
    
    def test_oligo_fasta_created(self):
        self.assertTrue(os.path.exists('oligo_seqs.fa'))
        
    def test_oligo_fasta_format_first_line(self):
        self.assertEqual(self.lines[0][0], '>')
        
    def test_oligo_fasta_format_last_line(self):
        self.assertTrue(bool(re.match('[ACGT]', self.lines[-1][0])))
        
    def test_oligo_fasta_line_number_is_even(self):
        self.assertEqual(len(self.lines) % 2, 0)
        
    def test_oligo_size(self):
        self.assertEqual(len(self.lines[1]), self.oligo)
        
    def test_first_oligo_within_range(self):
        first_line = self.lines[0].strip('>')
        parts = re.split('\W+', first_line)
        self.assertGreaterEqual(int(parts[1]), self.start)

if __name__ == '__main__':
    unittest.main()
    