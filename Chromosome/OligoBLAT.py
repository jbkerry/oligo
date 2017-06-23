#!usr/bin/python

from __future__ import division, print_function
import getopt,sys,re,datetime,pandas as pd

# Functions
    
def usage():
    print("usage: OligoBLAT.py -i <oligo numbers> -c <chromosome number>")
    
def GCcontent(x):
    gc_perc = (x.count('C') + x.count('G'))/len(x)
    return gc_perc

suffix = ""
chr_num = ""
BLAT=0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:c:h',)
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
        elif opt == '-i':
            suffix = arg
        elif opt == '-c':
            chr_num = arg
        else:
            usage()
            sys.exit(2)

chr_name = "chr"+chr_num
sys.stdout.write("Started running at: ")
TimeNow = datetime.datetime.now().time()
sys.stdout.write(str(TimeNow)+"\n")

# Main dictionaries

GC_dict = {}
AllOligos = {}
Deletions = {}
Matches = {}
Sequences = {}
Groups = {}

depthgauge = {}
Maffed = {}
used = {}
oligoCoverage = {}
OligoValue = {}
OligoSize = {}
StoredSeq = {}
dup_counts = {}
ProblemSequences = {}
qGapBases = {}

MafLine = 0
MafHeader = ""
MafSeq = ""

MafFile = "../"+chr_name+"_Oligos_"+suffix+".fa"

MafLines = [MafLineRead.rstrip('\n') for MafLineRead in open(MafFile)]
for ThisMafLine in MafLines:    
    MafLine=MafLine+1
    if MafLine==1:
        MafHeader = ThisMafLine
        MafHeader = re.sub('>','',MafHeader)
        
        Mchr, Mstart, Mstop, Mfragstart, Mfragstop, Mside = re.split("\W+",MafHeader)
        if Mchr not in depthgauge.keys():
            depthgauge[Mchr] = {}
        if Mchr not in oligoCoverage.keys():
            oligoCoverage[Mchr] = {}
        Mcounter = int(Mstart)
        Maffed[MafHeader]=0
        while Mcounter<=int(Mstop):
            if Mcounter not in depthgauge[Mchr].keys():
                depthgauge[Mchr].update({Mcounter: 0})
            if Mcounter not in oligoCoverage[Mchr].keys():
                oligoCoverage[Mchr].update({Mcounter: 1})
            else:
                oligoCoverage[Mchr][Mcounter]=oligoCoverage[Mchr][Mcounter]+1
            Mcounter=Mcounter+1
    elif MafLine==2:
        MafSeq = ThisMafLine
        StoredSeq[MafHeader]=MafSeq
        gc_perc = GCcontent(MafSeq)
        GC_dict[MafHeader] = gc_perc
        MafLine=0       

# Determine alignment matches and mismatches
BLATfile = './'+chr_name+'_'+suffix+'.psl'
linenumber = 0;
BLATLines = [BLATLine.rstrip('\n') for BLATLine in open(BLATfile)]
for ThisBLATLine in BLATLines[5:]:  

    #linenumber=linenumber+1
       
    match, mismatch, repmatch, ns, qgap, qgapbases, tgap, tgapbases, strand, query, qsize, qstart, qend, tname, tsize, tstart, tend, blockcount, blocksize, qtarts, tstarts = re.split("\s+",ThisBLATLine)
    if re.search('_',MafHeader)!=None:
        next
       
    percent = int(qsize) / 100      
    percent_match = int(match) / percent  
    if percent_match>=70:
        if query not in dup_counts.keys():
            dup_counts[query]=1
        else:
            dup_counts[query]=dup_counts[query]+1
    
    if query not in qGapBases.keys():
        qGapBases[query]=int(qgapbases)
    else:
        qGapBases[query]=qGapBases[query]+int(qgapbases)
    Oligochr, start, stop, fragstart, fragstop, side=re.split("\W+",query)            
        
    if query not in used.keys():
        counter = int(start)
        GoUntil = int(start)+int(qsize)
        while counter <= GoUntil:
            oligoCoverage[Oligochr][counter]=oligoCoverage[Oligochr][counter]+1;
            counter=counter+1 
       
    Maffed[query]=Maffed[query]+1
    if query not in used.keys():
        used[query]=1
    else:
        used[query]=used[query]+1
    OligoValue[query] = 0           
       
    ### This is where it tots up all blat hits for a given oligo
       
    loopstart = int(start) + int(qstart)
    counter = int(qstart)
       
    while counter<=int(qend):
        depthgauge[Oligochr][loopstart]=depthgauge[Oligochr][loopstart]+1
        counter=counter+1
        loopstart=loopstart+1         


for StoredID in Maffed.keys():   
    storedchr, storedstr, storedstp, storedFragstr, storedFragstp, storedSide=re.split("\W+",StoredID); 
    posCounter = int(storedstr)

    size = int(storedstp) - int(storedstr)
    OligoSize[StoredID]=size

    while posCounter <= int(storedstp):
        if posCounter in depthgauge[storedchr].keys():
            OligoValue[StoredID]=OligoValue[StoredID]+depthgauge[storedchr][posCounter]
            posCounter=posCounter+1
        else:
            print("Warning: entry missing from blat file")
    OligoValue[StoredID]=OligoValue[StoredID]-qGapBases[StoredID]

# Repeat Masker

SSRLength_dict = {}
SSRType_dict = {}
RM_file = "./"+chr_name+"_Oligos_"+suffix+".fa.out"
RMlines = [RMline.rstrip('\n') for RMline in open(RM_file)]
for ThisRMline in RMlines[3:]:
    parts = re.split("\s+",ThisRMline)
    Qname = parts[5]
    RepeatType = parts[10]
    Chr,Start,Stop,FragStart,FragEnd,Side = re.split("\W+",Qname)
    
    if len(Side)>1: # This is for oligos that have more than one repeat in them
        Side=Side[0]
        Qname = Chr+":"+Start+"-"+Stop+"-"+FragStart+"-"+FragEnd+"-"+Side
        
    Qstart = int(parts[6])
    Qstop = int(parts[7])
    SSRLength = (Qstop - Qstart)+1
    if Qname not in SSRLength_dict.keys():
        SSRLength_dict[Qname] = SSRLength # Store SSR length as this will be used to filter oligos
        SSRType_dict[Qname] = RepeatType
    else:
        if SSRLength>SSRLength_dict[Qname]:
            SSRLength_dict[Qname]=SSRLength
            SSRType_dict[Qname] = RepeatType
        
# Write text file with oligo info
Written = {}
OligoFile = "OligoInfo.txt"
#OligoFile = "OligoInfo_"+suffix+".txt"
TextOut = open(OligoFile,"w")
TextOut.write("Chr\tStart\tStop\tFragment Start\tFragment Stop\tSide of fragment\tSequence\tDensity score\tRepeat length\tRepeat Class\tGC%\n")
for ThisOligo in OligoValue.keys():
    density = OligoValue[ThisOligo]/OligoSize[ThisOligo]
    density_round = float("{0:.2f}".format(density))
    Write = 0
    Chr,Start,Stop,FragStart,FragEnd,Side = re.split("\W+",ThisOligo)
    OligoCoor = Chr+":"+Start+"-"+Stop
    if OligoCoor not in Written.keys():
        Write=1
    RepeatLength = 0
    RepeatType = "NA"
    if ThisOligo in SSRLength_dict.keys():
        RepeatLength=SSRLength_dict[ThisOligo]
        RepeatType=SSRType_dict[ThisOligo]
    if Write==1:   
        TextOut.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.2f}\t{8}\t{9}\t{10:.2f}\n".format(Chr,Start,Stop,FragStart,FragEnd,Side,StoredSeq[ThisOligo],density_round,RepeatLength,RepeatType,GC_dict[ThisOligo]))
        Written[OligoCoor] = 1
TextOut.close()

sys.stdout.write("Finished running at: ")
TimeNow = datetime.datetime.now().time()
sys.stdout.write(str(TimeNow)+"\n")
