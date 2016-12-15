<style>
u {
    text-decoration: underlne;
}
</style>

# Whole Chromosome Oligo Design

<u>Description</u><br>
This pipeline will generate all the unique oligos adjacent to a specified restriction site for a given chromosome. The user is required to specify the genome build, chromosome number, restricition enzyme and size of the oligos.
The pipeline provides the user with information about off-target binding, the presence of simple-sequence repeats and GC content for every oligo.

<u>Input</u><br>
The entire pipe can be run by simply supplying WholeChrPipe.sh with the variables <b>Genome</b>, <b>Chr</b>, <b>Enzyme</b> and <b>Oligo</b>.<br>
<b>Genome</b>: select from <b>hg18</b>, <b>hg19</b>, <b>mm9</b> or <b>mm10</b><br>
<b>Chr</b>: supply just the chromosome number or letter e.g. <b>7</b> or <b>X</b><br>
<b>Enzyme</b>: choose from <b>DpnII</b> (GATC), <b>NlaIII</b> (CATG) or <b>HindIII</b> (AAGCTT)<br>
<b>Oligo</b>: supply the number of bp for the required oligo length e.g. <b>70</b><br>

Example run for 70bp oligos adjacent to DpnII restriction sites on mouse mm9 chromosome 11:

<b>bash WholeChrPipe.sh Genome=mm9,Chr=11,Enzyme=DpnII,Oligo=70</b>

All supplied arguments are case sensitive.

<u>Output</u><br>
A file called <b>AllOligos_info.txt</b> will be generated in the directory created by the pipe for this run, named e.g. <b>mm9_chr11_DpnII_70bp/</b>. This will contain information about all of the oligos.

<u>Under the hood</u><br>
Below is a breakdown of the pipeline workflow.

<b>Workflow of WholeChrPipe.sh:</b>
<table>
    <tr>
        <th>Order</th>
        <th>Script</th>
        <th>Subprocesses</th>
        <th>Run in:</th>
    </tr>
    <tr>
        <td align="center">1</td>
        <td>CaptureOligos_WholeChromosome.py</td>
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
        <td>QsubShells.sh</td>
        <td>STAR; RepeatMasker; STAR_density_v1-10.py</td>
        <td>Directory generated from SplitFA.sh e.g mm9_chr11_DpnII_70bp/</td>
    </tr>
    <tr>
        <td align="center">5</td>
        <td>cat_OligoInfo.sh</td>
        <td>n/a</td>
        <td>Directory generated from SplitFA.sh e.g mm9_chr11_DpnII_70bp/</td>
    </tr>
</table>
