#!usr/bin/bash
module load python
module load bedtools

while getopts ":b:g:o:s:d:h" opt; do
  case $opt in
    b) Bed="$OPTARG"
    ;;
    g) Genome="$OPTARG"
    ;;
    o) Oligo="$OPTARG"
    ;;
    s) Step="$OPTARG"
    ;;
    d) MaxDist="$OPTARG"
    ;;
    h) Help=1
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
    exit 1
    ;;
  esac
done

if [ -z "$Bed" ] || [ -z "$Genome" ] || [ -z "$Oligo" ] || [ -z "$Step" ] || [ -z "$MaxDist" ] && [ -z "$Help" ]
then
    echo "ERROR: All arguments must be supplied"
    echo "Usage is: bash OT_Pipe.sh -b <bed file> -g <genome> -o <oligo size (bp)> -s <step size (bp)> -d <maximum distance (bp)>"
    echo "For more information type: bash OT_Pipe.sh -h"
    exit 1
elif [ -n "$Help" ]
then
    echo -e "\nOligo Design - CRISPR off-target\n"
    echo -e "-------------------------------------------------------\n"
    echo -e "This pipeline will generate oligos for performing Capture-C adjacent to predicted off-target cut sites of CRISPR. The user supplies a bed file containing the predicted off-target sites and oligos are generated in a step-wise manner walking away from the cut site in order to obtain the most efficient oligos within a required distance. The user can specify the oligo size, step size and maximum distance away from the cut site that the oligos are designed within.\n"
    
    echo -e "The pipeline can be run by supplying OT_Pipe.sh with the variables -b <bed file> -g <genome> -o <oligo size (bp)> -s <step size (bp)> and -d <maximum distance from cut site (bp)>."
    echo -e "-b supply the file name of a 4-column (Chr, Start, Stop, Name) bed file containing predicted off-target sites"
    echo -e "-g select from hg18, hg19, hg38, mm9 and mm10"
    echo -e "-o choose the size of the oligos (in bp) to be generated"
    echo -e "-s choose the step size (in bp) to specify the distance between adjacent oligos that are generated"
    echo -e "-d choose the maximum distance (in bp) away from the off-target site to design oligos\n"
    echo -e "-h print help and exit script\n"
    
    echo -e "Example run for 50bp oligos generated in a 10-bp stepwise manner no further than 200bp away from the off-target site, on either side (i.e. a maximum possible window of 400bp), for human hg19:\n"
    echo -e "bash OT_Pipe.sh -b OffTargetSites.bed -g hg19 -o 50 -s 10 -d 200\n"
    echo -e "All supplied arguments are case sensitive\n"
    exit 1
fi
if egrep -q '%|_' TSSCoors.bed; then
    echo "percent symbols (%) and underscores (_) are not allowed in gene names or coordinates"
    exit 1
fi
organism=""
species=""
if [ $Genome == "hg18" ] || [ $Genome == "hg19" ] || [ $Genome == "hg38" ]; then
    organism="Homo_sapiens"
    species="human"
elif [ $Genome == "mm9" ] || [ $Genome == "mm10" ]; then
    organism="Mus_musculus"
    species="mouse"
fi
echo "Generating possible oligo coordinates..."
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/CRISPR-OT/OT_OligoGen.py -b $Bed -g $Genome -o $Oligo -s $Step -d $MaxDist
echo "Extracting oligo sequences..."
fastaFromBed -fi /databank/igenomes/$organism/UCSC/$Genome/Sequence/WholeGenomeFasta/genome.fa -bed ./OToligoCoor.bed -name -fo OT_seqs.fa
echo "Running sequences through STAR..."
/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ./OT_seqs.fa --genomeDir /databank/igenomes/$organism/UCSC/$Genome/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix OT_Oligos_
echo "Checking for repeats..."
repeatmasker -noint -s -species $species ./OT_seqs.fa
echo "Selecting most efficient oligos..."
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/CRISPR-OT/OT_STAR.py -b $Bed
sort -k1,1 -k2,2n OT_OligoInfo.txt >AllOligos_Info.txt
sort -k1,1 -k2,2n OT_OligoInfo_filtered.txt >FilteredOligos_Info.txt
rm -f OT_OligoInfo.txt
rm -f OT_OligoInfo_filtered.txt
