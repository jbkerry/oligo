#!usr/bin/bash
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/CaptureOligos_WholeChromosome.py -g $Genome -c $Chr -e $Enzyme -o $Oligo
FastaName="Oligos_"$Genome"_chr"$Chr"_"$Enzyme"_"$Oligo"bp.fa"
DirName=${FastaName:7:$((${#FastaName}-10))}
lines=$(wc -l < ../$FastaName) # store number of lines
Max=$((((($lines/2)/20000)+1)*20000))
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/SplitFA.sh $FastaName
cd $DirName
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/MakeShells.py -g $Genome -c $Chr -u $Max
mkdir Logs
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/QsubShells.sh $Max