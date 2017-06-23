#!/usr/bin/env python

import getopt,sys,math

def usage():
    print("usage: OT_OligoGen.py -b <bed file> -g <genome build> -o <oligo size (bp)> -s <step size (bp)> -d <max distance from cut-site (bp)>")

bedfile = ""
genome = ""
oligo_size = 0
step_size = 0
max_distance = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'b:g:o:s:d:h',)
except getopt.GetoptError:
    usage()
    sys.exit(2)
    
if not opts:
    usage()
    sys.exit(2)
else:
    for opt, arg in opts:
        if opt == '-h':
            help_info()
            sys.exit(2)
        elif opt == '-b':
            bedfile = arg
        elif opt == '-g':
            genome = arg
        elif opt == '-o':
            oligo_size = arg
        elif opt == '-s':
            step_size = arg
        elif opt == '-d':
            max_distance = arg
        else:
            usage()
            sys.exit(2)
            
KillScript=0
if (genome!="hg18") & (genome!="hg19") & (genome!="hg38") & (genome!="mm9") & (genome!="mm10"):
    print("Genome not recognised, please choose from hg18, hg19, hg38, mm9 or mm10 (case sensitive)")
    KillScript=1

try:
    oligo_size = int(oligo_size)
    if oligo_size<1:
        print("Oligo size invalid, it must be an integer greater than 0")
        KillScript=1
except ValueError:
    print("Oligo size invalid, it must be an integer greater than 0")
    KillScript=1
    
try:
    step_size = int(step_size)
    if step_size<1:
        print("Step size invalid, it must be an integer greater than 0")
        KillScript=1
except ValueError:
    print("Step size invalid, it must be an integer greater than 0")
    KillScript=1
    
try:
    max_distance = int(max_distance)
    if max_distance<1:
        print("Maximum distance invalid, it must be an integer greater than 0")
        KillScript=1
except ValueError:
    print("Maximum distance invalid, it must be an integer greater than 0")
    KillScript=1
    
if KillScript==1:
    sys.exit(2)

OligoNumber = math.floor((max_distance-oligo_size)/step_size)
#OligoNumber = 15

#BedFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/GWAS/PipeTest/TSSCoors.bed"
OutFile = open("OToligoCoor.bed","w")
BedLines = [BedLine.rstrip('\n') for BedLine in open(bedfile)]
for ThisLine in BedLines:
    Chr, Start, Stop, Name = ThisLine.split('\t')
    Counter=1
    ProbeUpEnd = int(Start)-10
    ProbeUpStart = ProbeUpEnd-oligo_size
    ProbeDownStart = int(Stop)+10
    ProbeDownEnd = ProbeDownStart+oligo_size
    while Counter<=OligoNumber:
        #ProbeUpName = Chr+":"+ProbeUpStart+"-"+ProbeUpEnd+"%"+Name+"_Up"+Counter
        ProbeUpName = "{0}:{1}-{2}%{3}_Up{4}".format(Chr,ProbeUpStart,ProbeUpEnd,Name,Counter)
        ProbeDownName = "{0}:{1}-{2}%{3}_Down{4}".format(Chr,ProbeDownStart,ProbeDownEnd,Name,Counter)
        #ProbeDownName = Chr+":"+ProbeDownStart+"-"+ProbeDownEnd+"%"+Name+"_Down"+Counter
        #ProbeUpName =~ s/\r//g
        #ProbeDownName =~ s/\r//g
        OutFile.write("{0}\t{1}\t{2}\t{3}\n".format(Chr,ProbeUpStart,ProbeUpEnd,ProbeUpName))
        OutFile.write("{0}\t{1}\t{2}\t{3}\n".format(Chr,ProbeDownStart,ProbeDownEnd,ProbeDownName))
        ProbeUpStart-=step_size
        ProbeUpEnd-=step_size
        ProbeDownStart+=step_size
        ProbeDownEnd+=step_size
        Counter+=1
OutFile.close()