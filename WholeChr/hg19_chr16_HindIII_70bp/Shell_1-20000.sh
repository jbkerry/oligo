mkdir 1-20000
cd 1-20000
/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ../chr16_Oligos_1-20000.fa --genomeDir /databank/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr16_1-20000_
repeatmasker -noint -s -dir ./ -species human ../chr16_Oligos_1-20000.fa
python /t1-home/nuffmed/jkerry/Python/STAR_density_v1-10.py -i 1-20000 -c 16
