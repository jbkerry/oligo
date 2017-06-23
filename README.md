# Capture-C Oligo Design
This work consists of three modules for Capture-C and FISH oligo design.<br>
### __Regions__
This module will generate oligos adjacent to the nearest restriction sites flanking user-supplied coordinates for a given restriction enzyme. The user is required to specify a 4-column bed file (in the format 'Chr'\t'Start'\t'Stop'\t'e.g.Gene Name'), genome build, restriction enzyme,
oligo length and whether to run through BLAT or STAR. Functions in the 'tools' module can then be used to provide the user with information about off-targets binding and the presence of simple-sequence repeats within the oligo sequences. The full pipeline can be run from the command
line as shown below<br>
```
usage: regions.py [-h] -f FASTA -g GENOME -b BED [-o OLIGO] [-e ENZYME]<br>
                  [-s STAR_INDEX] [--blat]<br><br>

optional arguments:<br>
  -h, --help            show this help message and exit<br>
  -f FASTA, --fasta FASTA<br>
                        Path to reference genome fasta.<br>
  -g GENOME, --genome GENOME<br>
                        Genome build e.g. 'mm10' or 'hg38'.<br>
  -b BED, --bed BED     Path to bed file with capture viewpoint coordinates<br>
  -o OLIGO, --oligo OLIGO<br>
                        The size (in bp) of the oligo to design, default=70<br>
  -e ENZYME, --enzyme ENZYME<br>
                        Name of restriction enzyme, default=DpnII<br>
  -s STAR_INDEX, --star_index STAR_INDEX<br>
                        Path to STAR index directory. Omit this option if<br>
                        running with BLAT (--blat)<br>
  --blat                Detect off-targets using BLAT instead of STAR.<br><br>
  
e.g.<br><br>

python regions.py -f ~/mm9/genome.fa -g mm9 -b ~/files/coordinates.bed -o 70 -e DpnII -s ~/mm9/STAR<br><br>
Depends on tools.py

```<br><br>
### <a href="https://github.com/jbkerry/OligoDesign/tree/master/Regions">__Regions__</a>
Generates oligos for user-supplied coordinates across the entire genome<br><br>
### <a href="https://github.com/jbkerry/OligoDesign/tree/master/CRISPR-OT">__CRISPR-OT__</a>
Generates oligos for checking predicted off-target sites of CRISPR cutting<br>

