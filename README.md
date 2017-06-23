# Capture-C Oligo Design
This work consists of three modules for Capture-C and FISH oligo design.<br>
### __Regions__
This module will generate oligos adjacent to the nearest restriction sites flanking user-supplied coordinates for a given restriction enzyme. The user is required to specify a 4-column bed file (in the format 'Chr'\t'Start'\t'Stop'\t'e.g.Gene Name'), genome build, restriction enzyme,
oligo length and whether to run through BLAT or STAR. Functions in the 'tools' module can then be used to provide the user with information about off-targets binding and the presence of simple-sequence repeats within the oligo sequences. The full pipeline can be run from the command
line, as shown below.<br>
```
usage: regions.py [-h] -f FASTA -g GENOME -b BED [-o OLIGO] [-e ENZYME]
                  [-s STAR_INDEX] [--blat]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to reference genome fasta.
  -g GENOME, --genome GENOME
                        Genome build e.g. 'mm10' or 'hg38'.
  -b BED, --bed BED     Path to bed file with capture viewpoint coordinates
  -o OLIGO, --oligo OLIGO
                        The size (in bp) of the oligo to design, default=70
  -e ENZYME, --enzyme ENZYME
                        Name of restriction enzyme, default=DpnII
  -s STAR_INDEX, --star_index STAR_INDEX
                        Path to STAR index directory. Omit this option if
                        running with BLAT (--blat)
  --blat                Detect off-targets using BLAT instead of STAR.
  
e.g.

python regions.py -f ~/mm9/genome.fa -g mm9 -b ~/files/coordinates.bed -o 70 -e DpnII -s ~/mm9/STAR

Depends on tools.py
```

### <a href="https://github.com/jbkerry/OligoDesign/tree/master/Regions">__Regions__</a>
Generates oligos for user-supplied coordinates across the entire genome<br><br>
### <a href="https://github.com/jbkerry/OligoDesign/tree/master/CRISPR-OT">__CRISPR-OT__</a>
Generates oligos for checking predicted off-target sites of CRISPR cutting<br>

