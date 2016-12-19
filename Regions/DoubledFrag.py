#!/usr/bin/env python

#OligoFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/Oligos_nonDHS.txt"
OligoFile = "./Oligos_filtered.txt"
OligoLines = [OligoLine.rstrip('\n') for OligoLine in open(OligoFile)]
DoubleCounter = 0
FragCounter = 0
Dict={}
out_file = open("stats.txt","w")
for ThisOligoLine in OligoLines[1:]:
    parts = ThisOligoLine.split('\t')
    FragCoor = parts[0]+":"+parts[3]+"-"+parts[4]
    if FragCoor not in Dict.keys():
        Dict[FragCoor]=1
        FragCounter = FragCounter+1
    else:
        DoubleCounter=DoubleCounter+1
out_file.write("Following BLAT and RM filter there are "+str(FragCounter)+" fragments in total, "+str(DoubleCounter)+" of which have two oligos")