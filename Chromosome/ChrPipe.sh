#!usr/bin/bash
if [ -z "$Region" ]
then
    python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoGen.py -g $Genome -c $Chr -e $Enzyme -o $Oligo
    FastaName="Oligos_"$Genome"_chr"$Chr"_"$Enzyme"_"$Oligo"bp.fa"
else
    python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoGen.py -g $Genome -c $Chr -e $Enzyme -o $Oligo -r $Region
    FastaName="Oligos_"$Genome"_chr"$Chr"_"$Region"_"$Enzyme"_"$Oligo"bp.fa"
fi
if [ -z "$BLAT" ]
then
    BLAT=0
fi
TruncFasta=$(echo $FastaName | egrep -o "^[^.]+")
DirName=${TruncFasta:7}
##DirName=${FastaName:7:$((${#FastaName}-10))}
lines=$(wc -l < $FastaName) # store number of lines
MidMax="$(awk 'BEGIN { rounded = sprintf("%.0f", '$lines'/40000+0.4999999999); print rounded }')"
Max=$(($MidMax*20000))
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/SplitFA.sh $FastaName
cd $DirName
python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/MakeShells.py -g $Genome -c $Chr -u $Max -b $BLAT
mkdir Logs
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/RunShells.sh $Max
if [ $BLAT == 0 ]
then
    echo -e "Chr\tStart\tStop\tFragment Start\tFragment Stop\tSide of fragment\tSequence\tTotal number of alignments\tSTAR Density score\tRepeat length\tRepeat Class\tGC%" >AllOligos_info.txt
else
    echo -e "Chr\tStart\tStop\tFragment Start\tFragment Stop\tSide of fragment\tSequence\tBLAT Density score\tRepeat length\tRepeat Class\tGC%" >AllOligos_info.txt
fi
bash /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoMerge.sh $Max


##while getopts ":a:p:" opt; do
##  case $opt in
##    a) arg_1="$OPTARG"
##    ;;
##    p) p_out="$OPTARG"
##    ;;
##    \?) echo "Invalid option -$OPTARG" >&2
##    ;;
##  esac
##done
##
##printf "Argument p_out is %s\n" "$p_out"
##printf "Argument arg_1 is %s\n" "$arg_1"
##Then you can do
##
##$ ./my_script -p '/some/path' -a5
##Argument p_out is /some/path
##Argument arg_1 is 5