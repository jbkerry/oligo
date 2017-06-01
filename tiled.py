#!/usr/bin/env python

class Capture(object):
    '''Designs oligos for capture from adjacent restriction sites within a
    user-specified region
    
    Parameters
    ----------
    genome: hg18, hg19, hg38, mm9 or mm10
    chromosome: chromosome number/letter e.g. 7 or X
    enzyme: DpnII (GATC), NlaIII (CATG) or HindIII (AAGCTT)
    oligo: the length of the oligo to design (bp)
    blat: boolean, check off-targets using BLAT instead of STAR, default=False
    '''
    
    def __init__(self, genome, chromosome, enzyme, oligo, blat=False):
        self.genome = genome
        self.chromosome = chromosome
        self.enzyme = enzyme
        self.oligo = oligo
        self.blat = blat
    
    def generate_oligos():
        pass
