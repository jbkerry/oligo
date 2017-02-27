#!/usr/bin/env python

OligoSize = 50
StepSize = 10
OligoNumber = 15

BedFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/GWAS/PipeTest/TSSCoors.bed"
OutFile = open("OToligoCoor.bed","w")
BedLines = [BedLine.rstrip('\n') for BedLine in open(BedFile)]
for ThisLine in BedLines:
    Chr, Start, Stop, Name = ThisLine.split('\t')
    Counter=1
    ProbeUpEnd = int(Start)-10
    ProbeUpStart = ProbeUpEnd-OligoSize
    ProbeDownStart = int(Stop)+10
    ProbeDownEnd = ProbeDownStart+OligoSize
    while Counter<=OligoNumber:
        #ProbeUpName = Chr+":"+ProbeUpStart+"-"+ProbeUpEnd+"%"+Name+"_Up"+Counter
        ProbeUpName = "{0}:{1}-{2}%{3}_Up{4}".format(Chr,ProbeUpStart,ProbeUpEnd,Name,Counter)
        ProbeDownName = "{0}:{1}-{2}%{3}_Down{4}".format(Chr,ProbeDownStart,ProbeDownEnd,Name,Counter)
        #ProbeDownName = Chr+":"+ProbeDownStart+"-"+ProbeDownEnd+"%"+Name+"_Down"+Counter
        #ProbeUpName =~ s/\r//g
        #ProbeDownName =~ s/\r//g
        OutFile.write("{0}\t{1}\t{2}\t{3}\n".format(Chr,ProbeUpStart,ProbeUpEnd,ProbeUpName))
        OutFile.write("{0}\t{1}\t{2}\t{3}\n".format(Chr,ProbeDownStart,ProbeDownEnd,ProbeDownName))
        ProbeUpStart-=StepSize
        ProbeUpEnd-=StepSize
        ProbeDownStart+=StepSize
        ProbeDownEnd+=StepSize
        Counter+=1
OutFile.close()