#!usr/bin/bash
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/CaptureOligos_WholeChromosome.py -g mm9 -c 11 -e HindIII -o 70
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/SplitFA.sh Oligos_mm9_chr11_HindIII_70bp.fa
cd mm9_chr11_HindIII_70bp/
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/MakeShells.py -g mm9 -c 11 -u 80000
mkdir Logs
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/WholeChr/QsubShells.sh 80000