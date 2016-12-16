#!usr/bin/bash
if [ -z "$Region" ]
then
    python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoGen.py -g $Genome -c $Chr -e $Enzyme -o $Oligo
else
    python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoGen.py -g $Genome -c $Chr -e $Enzyme -o $Oligo -r $Region
fi
FastaName="Oligos_"$Genome"_chr"$Chr"_"$Enzyme"_"$Oligo"bp.fa"
DirName=${FastaName:7:$((${#FastaName}-10))}
lines=$(wc -l < $FastaName) # store number of lines
Max=$((((($lines/2)/20000)+1)*20000)) # This needs to be corrected for if the number of lines/2 equals a number exactly divisible by 20000
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/SplitFA.sh $FastaName
cd $DirName
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/MakeShells.py -g $Genome -c $Chr -u $Max
mkdir Logs
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/RunShells.sh $Max
echo -e "Chr\tStart\tStop\tFragment Start\tFragment Stop\tSide of fragment\tSequence\tTotal number of alignments\tDensity score\tRepeat length\tRepeat Class\tGC%" >AllOligos_info.txt
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoMerge.sh $Max