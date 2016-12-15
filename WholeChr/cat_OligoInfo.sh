#!usr/bin/bash

#touch AllOligos_info.txt
Counter=1
TopCounter=20000
Max=$1
Increment=20000
while [ $TopCounter -le $Max ]; do
    sed -e '1d' ./$Counter-$TopCounter/OligoInfo_$Counter-$TopCounter.txt >TempFile.txt
    cat AllOligos_info.txt >TempFile2.txt
    cat TempFile2.txt TempFile.txt >AllOligos_info.txt
    let Counter=Counter+Increment
    let TopCounter=TopCounter+Increment
done
rm -f TempFile*
