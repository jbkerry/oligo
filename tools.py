#!/usr/bin/env python

import os
import subprocess
import pysam
import re

path_list = [x.rstrip('\n') for x in open('config.txt')]
path_dict = dict(item.split(' = ') for item in path_list)

blat_param = '-stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999'
star_param = '--runThreadN 4 --genomeLoad NoSharedMemory ' \
             '--outFilterMultimapScoreRange 1000 --outFilterMultimapNmax ' \
             '100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 ' \
             '--seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax ' \
             '15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 ' \
             '--limitOutSAMoneReadBytes 300000 --outFileNamePrefix tiled_'

def check_off_target(fa, species, s_idx='', blat=False):
    '''Checks for repeat sequences in oligos generated from
    generate_oligos() using RepeatMasker and checks for off-target binding
    using either BLAT or STAR
    
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
    rm_out = open('rm_log.txt', 'w')
    subprocess.run(
        '{} -noint -s -species {} {}'.format(rm_path, species, 'oligo_seqs.fa'),
        shell = True,
        stdout = rm_out,
        stderr = rm_out,
    )
    rm_out.close()
    
    if blat:
        path = os.path.join(path_dict['BLAT_PATH'], 'blat')
        print("Checking off-target binding with BLAT...")
        blat_out = open('blat_log.txt', 'w')
        subprocess.run(
            '{} {} {} {} blat_out.psl'.format(path, blat_param,
                                              fa, 'oligo_seq.fa'),
            shell=True,
            stdout = blat_out,
            stderr = blat_out,
        )
        blat_out.close()
    else:
        path = os.path.join(path_dict['STAR_PATH'], 'STAR')
        print("Checking off-target binding with STAR...")
        star_out = open('star_log.txt', 'w')
        subprocess.run(
            '{} --readFilesIn {} --genomeDir {} {}'.format(path,
                                                           'oligo_seqs.fa',
                                                           s_idx,
                                                           star_param),
            shell = True,
            stdout = star_out,
            stderr = star_out,
        )
        star_out.close()
        
    return 'Off-target detection completed'

def _get_gc(x):
    gc_perc = (x.count('C') + x.count('G'))/len(x)
    gc_perc = float("{0:.2f}".format(gc_perc))
    return gc_perc

def get_density(sam='tiled_Aligned.out.sam',
                blat_file='blat_out.psl', blat=False):
    all_oligos = {}
    with open('oligo_seqs.fa') as f:
        for x in f:
            header = re.sub('>', '', x.rstrip('\n'))
            seq = next(f).rstrip('\n')
            all_oligos[header] = [seq, _get_gc(seq), 0, 0,
                                       0, 0, 0, 'NA']
    if blat:        
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
        sf = pysam.AlignmentFile(sam, 'r')
        for r in sf.fetch(until_eof=True):
            if all_oligos[r.query_name][2] == 0:
                 all_oligos[r.query_name][2] = r.get_tag('NH')
                    
            for block in r.cigartuples:
                if block[0]==0:
                    all_oligos[r.query_name][3]+=block[1]
                elif (block[0]==1) | (block[0]==2):
                    all_oligos[r.query_name][4]+=block[1]
        
    for o in all_oligos:
        score = all_oligos[o][3]-all_oligos[o][4]
        density = score/len(all_oligos[o][0])
        all_oligos[o][5] = float("{0:.2f}".format(density))      
    
    rm_msg = _get_repeats(all_oligos=all_oligos)
    write_msg = _write_file(all_oligos=all_oligos)
    return 'Density scores calculated'
    print(rm_msg)
    print(write_msg)
    
def _get_repeats(all_oligos):
    with open('oligo_seqs.fa.out') as f:
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
                    qname, dup = re.split('_', qname)
                
                qstart, qstop = map(int, (parts[6:8]))
                length = (qstop - qstart)+1
                if length>all_oligos[qname][6]:
                    all_oligos[qname][6:] = length, rep_type
            msg = 'Repeat scores calculated'
        else:
            msg = 'No repeats detected'
    
    return msg

def _write_file(all_oligos, file_name='oligo_info.txt'):
    with open(file_name, 'w') as f:
        f.write('Chr\tStart\tStop\tFragment Start\tFragment Stop\t' \
                'Side of fragment\tSequence\tTotal number of alignments\t' \
                'Density score\tRepeat length\tRepeat Class\tGC%\n')
        for key, idx in all_oligos.items():
            chr_name, start, stop, fragstart, fragstop, side = re.split(
                                                                '\W+', key)
            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        chr_name, start, stop, fragstart, fragstop, side,
                        idx[0], idx[2], idx[5], idx[6], idx[7], idx[1]))
    
    return 'Oligo data written to '+file_name