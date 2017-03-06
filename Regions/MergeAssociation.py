#!/usr/bin/env python

#GeneAssoc = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/GeneAssociations.txt"
GeneAssoc = "./GeneAssociations.txt"
#OligoFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/Oligos_nonDHS_BLATScores.txt"
#OligoFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/Oligos_filtered.txt"
OligoFile = "./Oligos_filtered.txt"
OligoAllFile = "./Oligos_all.txt"
#OligoFile = "./OligoInfo.txt"

Gene_dict = {}
GeneLines = [GeneLine.rstrip('\n') for GeneLine in open(GeneAssoc)]
for ThisGeneLine in GeneLines[1:]:
    Coor,Genes = ThisGeneLine.split('\t')
    Gene_dict[Coor] = Genes

#output = open("Oligos_nonDHS.txt","w")
output = open("AllOligos.txt","w")
filtoutput = open("FilteredOligos.txt","w")
OligoLines = [OligoLine.rstrip('\n') for OligoLine in open(OligoFile)]
#output.write("Chr\tOligo Start\tOligo Stop\tFragment Start\tFragment Stop\tSide of Fragment\tDensity\tRepeat Length\tSequence\tAssociated genes\n")
output.write(OligoLines[0]+"\tAssociated genes\n")
for ThisOligoLine in OligoLines[1:]:
    parts = ThisOligoLine.split('\t')
    FragCoor = parts[0]+":"+parts[3]+"-"+parts[4]
    output.write(ThisOligoLine+"\t"+Gene_dict[FragCoor]+"\n")
output.close()
filtoutput.close()
