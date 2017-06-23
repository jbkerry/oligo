# Capture-C Oligo Design
This work consists of three modules for Capture-C and FISH oligo design.<br>
### __regions__
This module will generate oligos adjacent to the nearest restriction sites flanking user-supplied coordinates for a given restriction enzyme. The user is required to specify a 4-column bed file (in the format 'Chr'\t'Start'\t'Stop'\t'e.g.Gene Name'), genome build, restriction enzyme,
oligo length and whether to run through BLAT or STAR. Functions in the 'tools' module can then be used to provide the user with information about off-target binding and the presence of simple-sequence repeats within the oligo sequences. The full pipeline can be run from the command
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

### __tiled__
This module will generate oligos adjacent to all instances of a specified restriction enzyme recognition sequence for an entire chromosome or region within the chromosome. The user is required to specify the genome build, chromosome number, restriction enzyme, size of the oligos and the region of the chromosome,
if applicable. This can also be run in FISH mode which is restriction enzyme-independent and designs oligos adjacent to each other based on a user-specified step-size. Functions in the 'tools' module can then be used to provide the user with information about off-target binding
and the presence of simple-sequence repeats within the oligo sequences. The full pipeline can be run from the command line, as shown below.<br>
```
usage: tiled.py [-h] [--fish] -f FASTA -g GENOME -c CHR [-e ENZYME]
                [-t STEP_SIZE] [-o OLIGO] [-r REGION] [-s STAR_INDEX] [--blat]

optional arguments:
  -h, --help            show this help message and exit
  --fish                Run in FISH mode (restriciton enzyme independent).
  -f FASTA, --fasta FASTA
                        Path to reference genome fasta.
  -g GENOME, --genome GENOME
                        Genome build e.g. 'mm10' or 'hg38'.
  -c CHR, --chr CHR     Chromosome number/letter on which to design the
                        oligos.
  -e ENZYME, --enzyme ENZYME
                        Name of restriction enzyme, default=DpnII. Omit this
                        option if running in FISH mode (--fish)
  -t STEP_SIZE, --step_size STEP_SIZE
                        Step size (in bp) between adjacent oligos when running
                        in FISH mode (--fish), default=70. Omit this option if
                        you are not using the --fish flag
  -o OLIGO, --oligo OLIGO
                        The size (in bp) of the oligo to design, default=70
  -r REGION, --region REGION
                        The region in which to design the oligos; must be in
                        the format 'start-stop' e.g. '10000-20000'. Omit this
                        option to design oligos across the entire chromosome.
  -s STAR_INDEX, --star_index STAR_INDEX
                        Path to STAR index directory. Omit this option if
                        running with BLAT (--blat)
  --blat                Detect off-targets using BLAT instead of STAR.
  
e.g.

python tiled.py -f ~/hg19/genome.fa -g hg19 -c 10 -e DpnII -o 70 -r 10520000-10620000 -s ~/hg19/STAR

Depends on tools.py
```

### __off_target__
This module will generate oligos for performing a capture pull-down adjacent to predicted off-target cut sites of CRISPR guide RNA. The user supplies a bed file containing the predicted off-target sites and oligos are generated in a step-wise manner walking away from the cut site.
Functions in the 'tools' module can then be used in order to determine the most efficient oligos within a required distance. The user can specify the oligo size, step size and maximum distance away from the cut site that the oligos are designed within.
The full pipeline can be run from the command line, as shown below.<br>
```
usage: off_target.py [-h] -f FASTA -g GENOME -b BED [-t STEP_SIZE] [-o OLIGO]
                     [-m MAX_DIST] [-s STAR_INDEX] [--blat]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Path to reference genome fasta.
  -g GENOME, --genome GENOME
                        Genome build e.g. 'mm10' or 'hg38'.
  -b BED, --bed BED     Path to bed file with off-target coordinates
  -t STEP_SIZE, --step_size STEP_SIZE
                        Step size (in bp) between adjacent oligos, default=10
  -o OLIGO, --oligo OLIGO
                        The size (in bp) of the oligo to design, default=70
  -m MAX_DIST, --max_dist MAX_DIST
                        The maximum distance away from the off-target site to
                        design oligos to, default=200
  -s STAR_INDEX, --star_index STAR_INDEX
                        Path to STAR index directory. Omit this option if
                        running with BLAT (--blat)
  --blat                Detect off-targets using BLAT instead of STAR.
  
e.g.

python off_target.py -f ~/hg38/genome.fa -g hg38 --b ~/files/offtargets.bed -o 70 -t 10 -m 200 --blat

Depends on tools.py
```

