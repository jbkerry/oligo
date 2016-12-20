#!/usr/bin/env python

import getopt,sys

def usage():
    print("usage: MakeShells.py -g <genome build> -c <chromosome number> -u <upper limit> -b <use BLAT? 1 = BLAT, 0 = STAR>")
    
genome = ""
chromosome = ""
upperLimit = 0
BLAT = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'g:c:u:b:h',)
except getopt.GetoptError:
    usage()
    sys.exit(2)
    
if not opts:
    usage()
    sys.exit(2)
else:
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(2)
        elif opt == '-g':
            genome = arg
        elif opt == '-c':
            chromosome = arg
        elif opt == '-u':
            upperLimit = arg
        elif opt == '-b':
            BLAT = arg
        else:
            usage()
            sys.exit(2)

# Genome
KillScript = 0
MaxChromosome = 19
Species = "mouse"
Organism = ""
if (genome=="hg18") | (genome=="hg19"):
    Organism="Homo_sapiens"
    Species="human"
    MaxChromosome = 22
elif (genome=="mm9") | (genome=="mm10"):
    Organism="Mus_musculus"
    Species="mouse"
    MaxChromosome = 19
else:
    print("Genome not recognised, please choose from hg18, hg19, mm9 or mm10 (case sensitive)")
    KillScript=1
    
# Chromosome
try:
    chr_value = int(chromosome)
    if (chr_value<1) | (chr_value>MaxChromosome):
        print("Selected chromosome number does not exist for the specified genome build")
        KillScript=1
except ValueError:
    if (chromosome!="X") & (chromosome!="Y") & (chromosome!="M"):
        print("Selected chromosome number does not exist for the specified genome build (case sensitive for X, Y and M)")
        KillScript=1
    
if KillScript==1:
    sys.exit(2)
            
CounterLow = 1
Counterhigh = 20000
while Counterhigh<=int(upperLimit):
    Suffix = str(CounterLow)+"-"+str(Counterhigh)
    out_file = open("Shell_"+Suffix+".sh","w")
    out_file.write("mkdir "+Suffix+"\n")
    out_file.write("cd "+Suffix+"\n")
    if BLAT==0:
        out_file.write("/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ../chr"+chromosome+"_Oligos_"+Suffix+".fa --genomeDir /databank/igenomes/"+Organism+"/UCSC/"+genome+"/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr"+chromosome+"_"+Suffix+"_\n")
    elif BLAT==1:
        out_file.write("module load blat")
        out_file.write("blat -stepSize=5 -minScore=10 -minIdentity=0 -repMatch=999999 /databank/igenomes/$organism/UCSC/$genome/Sequence/WholeGenomeFasta/genome.fa ../chr"+chromosome+"_Oligos_"+Suffix+".fa chr"+chromosome+"_"+Suffix+".psl")
    out_file.write("repeatmasker -noint -s -dir ./ -species "+Species+" ../chr"+chromosome+"_Oligos_"+Suffix+".fa\n")
    if BLAT==0:
        out_file.write("python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoSTAR.py -i "+Suffix+" -c "+chromosome+"\n")
    elif BLAT==1:
        out_file.write("python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoBLAT.py -i "+Suffix+" -c "+chromosome+"\n")
    CounterLow = CounterLow+20000
    Counterhigh = Counterhigh+20000
out_file.close()