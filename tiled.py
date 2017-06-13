#!/usr/bin/env python

import re
import os
import subprocess
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

org_dict = {'mm9': 'Mus_musculus',
            'mm10': 'Mus_musculus',
            'hg18': 'Homo_sapiens',
            'hg19': 'Homo_sapiens',
            'hg38': 'Homo_sapiens'}

rs_dict = {'DpnII': 'GATC',
           'NlaIII': 'CATG',
           'HindIII': 'AAGCTT'}

blat_param = '-stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999'
star_param = '--runThreadN 4 --genomeLoad NoSharedMemory ' \
             '--outFilterMultimapScoreRange 1000 --outFilterMultimapNmax ' \
             '100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 ' \
             '--seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax ' \
             '15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 ' \
             '--limitOutSAMoneReadBytes 300000 --outFileNamePrefix tiled_'

class Capture(object):
    '''Designs oligos for capture from adjacent restriction sites within a
    user-specified region
    
    Parameters
    ----------
    fa: path to reference genome/chromsome fasta
    blat: boolean, check off-targets using BLAT instead of STAR (not
    recommended for large designs), default=False
    '''
    
    def __init__(self, fa, blat=False):
        self.fa = fa
        self.blat = blat
    
    def generate_oligos(self, chromosome, enzyme='DpnII', oligo=70, region=''):
        '''Generates fasta file containing the oligos for a specific
        restriction enzyme
        
        Parameters
        ----------
        chromosome: chromosome number/letter e.g. 7 or X
        enzyme: DpnII (GATC), NlaIII (CATG) or HindIII (AAGCTT), default=DpnII
        oligo: the length of the oligo to design (bp), default=70
        region: the region of the chromosome to design oligos, must be in the
            format start-stop, e.g. 10000-20000. Omit this option to design
            oligos over the entire chromosome
            
        Output
        ------
        oligo_seqs.fa: a FASTA file containing sequences of all the oligos
        
        '''
        
        chr_name = 'chr'+str(chromosome)
        
        print('Loading reference fasta file...')
        seq_dict = SeqIO.to_dict(SeqIO.parse(self.fa, 'fasta'))
        seq = seq_dict[chr_name].seq.upper()
        print('Fasta file loaded. Generating oligos...')
            
        if not region:
            start = 0; stop = len(seq)
        else:
            start, stop = list(map(int, region.split('-')))
        
        p = re.compile(rs_dict[enzyme])
        pos_list = []
        for m in p.finditer(str(seq[start:stop])):
            pos_list.append(m.start()+start)
        
        cut_size = len(rs_dict[enzyme])
        fa_w = open('oligo_seqs.fa', 'w')    
        for i in range(len(pos_list)-1):
            j = i + 1
            frag_len = pos_list[j]-pos_list[i]+cut_size
            if (frag_len>=oligo):
                l_start = pos_list[i]
                l_stop = l_start+oligo
                l_seq = seq[l_start:l_stop]
                
                r_stop = pos_list[j]+cut_size
                r_start = r_stop-oligo
                r_seq = seq[r_start:r_stop]
                
                fa_w.write('>{0}:{1}-{2}-{1}-{3}-L\n{4}\n'.format(chr_name,
                                                                  l_start,
                                                                  l_stop,
                                                                  r_stop,
                                                                  l_seq))
                if frag_len>oligo:
                    fa_w.write('>{0}:{1}-{2}-{3}-{2}-R\n{4}\n'.format(chr_name,
                                                                      r_start,
                                                                      r_stop,
                                                                      l_start,
                                                                      r_seq))
        fa_w.close()
        
        return "Wrote oligos to oligo_seqs.fa"
    
    def _split_fa():
        f = open('oligo_seqs.fa')
        file_content = f.readlines()
        split = 40000
        start = 1
        stop = 20000
        for lines in range(0, len(file_content), split):
            output_data = file_content[lines:lines+split]
            with open('{0}-{1}.fa'.format(start, stop),'w') as f_out:
                f_out.write(''.join(output_data))
            start+=20000; stop+=20000
            
        return 'Split files per 20,000 oligos'
    
    def check_off_target(species, o_fa, s_idx=''):
        '''Checks for repeat sequences in oligos generated from generate_oligos()
        using RepeatMasker and checks for off-target binding using either BLAT
        or STAR
        
        Parameters
        ----------
        species: the general name of the species e.g. human or mouse. From the
            RepeatMasker help page:
            
            ############
            
            The species name must be a valid NCBI Taxonomy Database species
            name and be contained in the RepeatMasker repeat database.
            Some examples are:
    
              human
              mouse
              rattus
              "ciona savignyi"
              arabidopsis
    
            Other commonly used species:
    
            mammal, carnivore, rodentia, rat, cow, pig, cat, dog, chicken,
            fugu, danio, "ciona intestinalis" drosophila, anopheles, elegans,
            diatoaea, artiodactyl, arabidopsis, rice, wheat, and maize
            
            ############
            
        o_fa: the oligo fasta file generated from generate_oligos()
        s_idx: the directory containing the STAR index for this genome (not
            required if blat=True)
            
        Output
        ------
        oligo_seq.fa.out: RepeatMasker output file
        tiled_Aligned.out.sam (STAR): alignment file from STAR
        blat_out.psl (BLAT): off-target binding file from BLAT
        *a number of other files are produced from these programs but are not
            required in this pipeline. They are not deleted in case the user is
            interested in their contents.
            
        '''
        
        rm_path = os.path.join(path_dict['RM_PATH'], 'RepeatMasker')
        print('Checking for repeat sequences in oligos...')
        subprocess.run('{} -noint -s -species {} {}'.format(rm_path, species,
                                                            o_fa), shell=True)
        
        if self.blat:
            path = os.path.join(path_dict['BLAT_PATH'], 'blat')
            print("Checking off-target binding with BLAT...")
            subprocess.run('{} {} {} {} blat_out.psl'
                           .format(path, blat_param, self.fa, o_fa),
                           shell=True)
        else:
            path = os.path.join(path_dict['STAR_PATH'], 'STAR')
            print("Checking off-target binding with STAR...")
            subprocess.run('{} --readFilesIn {} --genomeDir {} {}'
                           .format(path, o_fa, s_idx, star_param), shell=True)
            
        return 'Off-target detection completed'
    
    def _get_gc(x):
        gc_perc = (x.count('C') + x.count('G'))/len(x)
        gc_perc = float("{0:.2f}".format(gc_perc))
        return gc_perc
    
    def get_density(sam='tiled_Aligned.out.sam', fa='oligo_seqs.fa',
                    blat_file='blat_out.psl'):
        # seq, gc, nh, matches, gaps, density
        all_oligos = {}
        if self.blat:
            with open(fa) as f:
                for x in f:
                    header = re.sub('>', '', x.rstrip('\n'))
                    seq = next(f).rstrip('\n')
                    all_oligos[header] = [seq, _get_gc(seq), 0, 0, 0, 0]
                    
            with open(blat_file) as f:
                for _ in range(5):
                    next(f)
                for x in f:
                    parts = re.split("\s+", x.rstrip('\n'))
                    query = parts[9]
                    qgapbases, qstart, qend = map(int, (parts[5], parts[11], 
                                                        parts[12]))
                    all_oligos[query][2]+=1
                    all_oligos[query][3]+=(int(qend)-int(qstart))+1
                    all_oligos[query][4]+=int(qgapbases)  
        else:
            all_oligos = {}
            sf = pysam.AlignmentFile(sam, 'r')
            for r in sf.fetch(until_eof=True):
                if r.query_name not in all_oligos:
                    seq = r.query_sequence
                    if r.is_reverse:
                        seq = str(Seq(seq).reverse_complement())
                    all_oligos[r.query_name] = [seq, _get_gc(seq), r.get_tag('NH'),
                                                0, 0, 0]
                        
                for block in r.cigartuples:
                    if block[0]==0:
                        all_oligos[r.query_name][3]+=block[1]
                    elif (block[0]==1) | (block[0]==2):
                        all_oligos[r.query_name][4]+=block[1]
            
        for o in all_oligos:
            score = all_oligos[o][3]-all_oligos[o][4]
            density = score/len(all_oligos[o][0])
            all_oligos[o][5] = float("{0:.2f}".format(density))      
        
        return all_oligos   
        
    def get_repeats():
        rm_dict = {}
        rm_file = "./oligo_seqs.fa.out"
        rm_lines = [x.rstrip('\n') for x in open(rm_file)]
        with open(rm_file) as f:
            for _ in range(3):
                next(f)
            for x in f:
                parts = re.split("\s+", x)
                qname = parts[5]
                rep_type = parts[10]
                chr_name, start, stop, fragstart, fragend, side = re.split(
                                                                "\W+", qname)
                if len(side)>1:
                    side = side[0]
                    qname = '{}:{}-{}-{}-{}-{}'.format(chr_name, start, stop,
                                                   fragstart, fragend, side)
                qstart, qstop = map(int, (parts[6:8]))
                length = (qstop - qstart)+1
                str_length = rm_dict.get(qname, [0, ''])[0]
                if length>str_length:
                    rm_dict[qname] = [length, rep_type]
        
        return rm_dict
    