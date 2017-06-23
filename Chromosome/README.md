# Oligo Design - Chromosome
`Description`<br>
This pipeline will generate all the oligos adjacent to a specified restriction site for an entire chromosome or region within the chromosome. The user is required to specify the genome build, chromosome number, restriction enzyme, size of the oligos and the region of the chromosome, if applicable.
The pipeline provides the user with information about off-target binding, the presence of simple-sequence repeats and GC content for every oligo.

`Input`<br>
The pipeline can be run by supplying ChrPipe.sh with the variables <b>Genome</b>, <b>Chr</b>, <b>Enzyme</b>, <b>Oligo</b>, <b>Region</b> and <b>BLAT</b>.<br>
<b>Genome</b>: select from <b>hg18</b>, <b>hg19</b>, <b>mm9</b> or <b>mm10</b><br>
<b>Chr</b>: supply just the chromosome number or letter e.g. <b>7</b> or <b>X</b><br>
<b>Enzyme</b>: choose from <b>DpnII</b> (GATC), <b>NlaIII</b> (CATG) or <b>HindIII</b> (AAGCTT)<br>
<b>Oligo</b>: supply the number of bp for the required oligo length e.g. <b>70</b><br>
<b>Region</b> (optional): supply coordinates to only generate oligos within a specific region of the chromosome (must be in the format Start-Stop) e.g. <b>3050000-5000000</b>. Omit this option to run the pipe on the entire chromosome.<br>
<b>BLAT</b> (optional): choose whether to test for off-target binding using either STAR or BLAT. Select <b>0</b> (to run STAR) or <b>1</b> (to run BLAT). If this option is omitted, the pipeline will run STAR by default. BLAT is not recommended for large designs, particularly on the human genome, because it becomes incredibly slow.<br>

Example run for 70bp oligos adjacent to DpnII restriction sites on mouse mm9 chromosome 11 (entire chromosome), checking off-target binding using BLAT:

<b>bash ChrPipe.sh Genome=mm9,Chr=11,Enzyme=DpnII,Oligo=70,BLAT=1</b>

Example run for 50bp oligos adjacent to HindIII restriction sites in the 10500000-12000000 region of human hg19 chromosome 5, checking off-target binding using STAR:

<b>bash ChrPipe.sh Genome=hg19,Chr=5,Enzyme=HindIII,Oligo=50,Region=105000000-12000000</b>

All supplied arguments are case sensitive.

`Output`<br>
A file called <b>AllOligos_info.txt</b> will be generated in the directory created by the pipe for this run, named e.g. <b>mm9_chr11_DpnII_70bp/</b>. This will contain information about all of the oligos.

`Under the hood`<br>
Below is a breakdown of the pipeline workflow.

<b>Workflow of ChrPipe.sh:</b>
<table>
    <tr>
        <th>Order</th>
        <th>Script</th>
        <th>Subprocesses</th>
        <th>Run in:</th>
    </tr>
    <tr>
        <td align="center">1</td>
        <td>OligoGen.py</td>
        <td>BioPython: SeqIO</td>
        <td>Any</td>
    </tr>
    <tr>
        <td align="center">2</td>
        <td>SplitFA.sh</td>
        <td>n/a</td>
        <td>Directory with generated fasta file</td>
    </tr>
    <tr>
        <td align="center">3</td>
        <td>MakeShells.py</td>
        <td>n/a</td>
        <td>Directory generated from SplitFA.sh e.g mm9_chr11_DpnII_70bp/
    </tr>
    <tr>
        <td align="center">4</td>
        <td>RunShells.sh</td>
        <td>STAR/BLAT; RepeatMasker; OligoSTAR.py/OligoBLAT.py</td>
        <td>Directory generated from SplitFA.sh e.g mm9_chr11_DpnII_70bp/</td>
    </tr>
    <tr>
        <td align="center">5</td>
        <td>OligoMerge.sh</td>
        <td>n/a</td>
        <td>Directory generated from SplitFA.sh e.g mm9_chr11_DpnII_70bp/</td>
    </tr>
</table>
</body>
</html>
