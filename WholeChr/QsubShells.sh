#!usr/bin/bash
#/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_1-20000.fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad LoadAndKeep --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr1_1-20000_
#/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_20001-40000.fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad LoadAndKeep --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr1_20001-40000_

#/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_1-100.fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad LoadAndKeep --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr1_1-100_
#/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_101-200.fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad LoadAndKeep --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr1_101-200_

#repeatmasker -noint -s -species mouse /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_1-20000.fa
#python STAR_depthGauge_v1-10.py -i 1-100
#qsub -cwd -o 1-20000_dg.out -e 1-20000_dg.err -N DG1 <./ShellScripts/Shell1.sh
#qsub -cwd -o 1-40000_dg.out -e 1-40000_dg.err -N DG2 <./ShellScripts/Shell2.sh

#qsub -cwd -o ./Logs/1-100_dg.out -e ./Logs/1-100_dg.err -N DG1 <./ShellScripts/Shell1.sh
#qsub -cwd -o ./Logs/101-200_dg.out -e ./Logs/101-200_dg.err -N DG2 <./ShellScripts/Shell2.sh

#qsub -cwd -o ./Logs/1-20000.out -e ./Logs/1-20000.err -N DG1 <./Shell_1-20000.sh
#qsub -cwd -o ./Logs/20001-40000.out -e ./Logs/20001-40000.err -N DG2 <./Shell_20001-40000.sh

Counter=600001
TopCounter=620000
Max=740000
Increment=20000
while [ $TopCounter -le $Max ]; do
    qsub -cwd -o ./Logs/$Counter-$TopCounter""_density.out -e ./Logs/$Counter-$TopCounter""_density.err -N Den$Counter <./density_$Counter-$TopCounter.sh
    let Counter=Counter+Increment
    let TopCounter=TopCounter+Increment
done