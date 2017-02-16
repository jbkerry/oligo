# Oligo Design - Regions
`Description`<br>
This pipeline will generate oligos adjacent to the nearest restriction sites flanking user-supplied coordinates for a given restriction enzyme. The user is required to specify a 4-column bed file (in the format 'Chr'\t'Start'\t'Stop'\t'e.g.Gene Name'), genome build, restriction enzyme,
oligo length and whether to run through BLAT or STAR. The pipeline provides the user with information about off-targets binding and the presence of simple-sequence repeats within the oligo sequences.

`Input`<br>
The pipeline can be by supplying RegionsPipe.sh with unnamed variables in the following order: <b>bed file</b>, <b>genome</b>, <b>enzyme</b>, <b>oligo size</b> and <b>BLAT/STAR</b><br>
<b>Bed file</b>: supply the name of a 4-column bed file in the format 'Chr'\t'Start'\t'Stop'\t'Info'
<b>Genome</b>: select from <b>hg18</b>, <b>hg19</b>, <b>mm9</b> or <b>mm10</b><br>
<b>Enzyme</b>: choose from <b>DpnII</b> (GATC), <b>NlaIII</b> (CATG) or <b>HindIII</b> (AAGCTT)<br>
<b>Oligo</b>: supply the number of bp for the required oligo length e.g. <b>70</b><br>
<b>BLAT/STAR</b>: supply <b>0</b> (for BLAT) or <b>1</b> (for STAR)<br><br>

Example run for 70bp oligos adjacent to DpnII restriction sites adjacent to coordinates in MyFile.bed for mm9, checking off-target effects with STAR:<br><br>
<b>bash RegionsPipe.sh MyFile.bed mm9 DpnII 70 1</b><br>

All supplied arguments are case-sensitive