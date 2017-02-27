module load python
module load bedtools
echo "Generating possible possible oligo coordinates..."
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/CRISPR-OT/OT_OligoGen.py
echo "Extracting oligo sequences..."
fastaFromBed -fi /databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed ./OToligoCoor.bed -name -fo OT_seqs.fa
echo "Running sequences through STAR..."
/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ./OT_seqs.fa --genomeDir /databank/igenomes/$organism/UCSC/$genome/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix OT_Oligos_
