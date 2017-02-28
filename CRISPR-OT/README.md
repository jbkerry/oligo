# Oligo Design - CRISPR off-target
`Description`<br>
This pipeline will generate oligos for performing Capture-C adjacent to predicted off-target cut sites of CRISPR. The user supplies a bed file containing the predicted off-target sites and oligos are generated in a step-wise manner walking away from the cut site in order to obtain
the most efficient oligos within a required distance. The user can specify the oligo size, step size and maximum distance away from the cut site that the oligos are designed within.

`Input`
The pipeline can be run by supplying OT_Pipe.sh with the variables -b [bed file] -g \<genome\> -o <oligo size (bp)> -s <step size (bp)> and -d <maximum distance from cut site (bp)>.<br>
<b>-b</b> supply the file name of a 4-column (Chr, Start, Stop, Name) bed file containing predicted off-target sites<br>
<b>-g</b> select from <b>hg18</b>, <b>hg19</b>, <b>hg38</b>, <b>mm9</b> and <b>mm10</b><br>
<b>-o</b> choose the size of the oligos (in bp) to be generated<br>
<b>-s</b> choose the step size (in bp) to specify the distance between adjacent oligos that are generated<br>
<b>-d</b> choose the maximum distance (in bp) away from the off-target site to design oligos<br><br>

Example run for 50bp oligos generated in a 10-bp stepwise manner no further than 200bp away from the off-target site, on either side (i.e. a maximum possible window of 400bp), for human hg19:<br><br>
<b>bash OT_Pipe.sh -b OffTargetSites.bed -g hg19 -o 50 -s 10 -d 200</b><br><br>
All supplied arguments are case sensitive

