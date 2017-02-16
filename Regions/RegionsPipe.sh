#!usr/bin/bash
module load blat
bedfile=$1
genome=$2
enzyme=$3
oligo=$4
STARvar=$5
organism=""
species=""
if [ $genome == "hg18" ] || [ $genome == "hg19" ]; then
    organism="Homo_sapiens"
    species="human"
elif [ $genome == "mm9" ] || [ $genome == "mm10" ]; then
    organism="Mus_musculus"
    species="mouse"
fi
DirName=$genome"_"$enzyme"_"$oligo"bp"

mkdir $DirName
cd $DirName

echo "Generating fragments..."
python ../FragExtract.py -b $bedfile -g $genome -e $enzyme -o $oligo

if [ $STARvar == 0 ]
then
    echo "Running through BLAT..."
    blat -stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999 /databank/igenomes/$organism/UCSC/$genome/Sequence/WholeGenomeFasta/genome.fa ./GeneratedOligos.fa GeneratedOligos.psl
elif [ $STARvar == 1 ]
then
    echo "Running through STAR..."
    /package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ./GeneratedOligos.fa --genomeDir /databank/igenomes/$organism/UCSC/$genome/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix GeneratedOligos_
else
    echo "Please choose 0 or 1 for running BLAT or STAR"
    exit 1
fi

echo "Repeat masker..."
repeatmasker -noint -s -species $species ./GeneratedOligos.fa

if [ $STARvar == 0 ]
then
    python ../DepthGauge.py
    python ../MergeAssociation.py
    python ../DoubledFrag.py
elif [ $STARvar == 1 ]
then
    python ../OligoSTAR.py
fi
echo "All done. Check stats.txt to see number of successful fragments and oligos"
