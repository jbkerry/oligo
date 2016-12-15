#!/usr/bin/bash
counter=1
OligoStart=1
topCounter=40000
increment=40000
OligoStop=20000
OligoIncrement=20000
#loopLimit=915661
FastaName=$1
ChrName=$(echo $FastaName | grep -o "chr[0-9]*[A-Z]*")
DirName=${FastaName:7:$((${#FastaName}-10))}
mkdir $DirName
cd $DirName
lines=$(wc -l < ../$FastaName) # store number of lines
loopLimit=$(($lines+1))
FinalTop=$((($lines/$increment+1)*$increment)) # This needs to be corrected for if the number of lines equals a number exactly divisible by 40000

## Set up a loop to generate fasta files for every 20,000 oligos (40,000 lines)

while [ $topCounter -lt $loopLimit ]; do
    awk 'FNR>='$counter' && FNR<='$topCounter ../$FastaName >$ChrName""_Oligos_$OligoStart-$OligoStop.fa
    let counter=counter+increment
    let topCounter=topCounter+increment
    let OligoStart=OligoStart+OligoIncrement
    let OligoStop=OligoStop+OligoIncrement
    if [ $topCounter -eq $FinalTop ]
    then
        let topCounter=$lines
    fi
done

