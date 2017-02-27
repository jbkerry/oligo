module load python
module load bedtools
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
#/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ./OT_seqs.fa --genomeDir /databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix OT_Oligos_
