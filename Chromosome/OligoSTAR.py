#!usr/bin/python

from __future__ import division, print_function
import getopt,sys,re,datetime,pandas as pd

# Functions
    
def usage():
    print("usage: OligoSTAR.py -i <oligo numbers> -c <chromosome number>")
    
def GCcontent(x):
    gc_perc = (x.count('C') + x.count('G'))/len(x)
    return gc_perc

suffix = ""
chr_num = ""

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

# Determine alignment matches and mismatches
SAMfile = './'+chr_name+'_'+suffix+'_Aligned.out.sam'
SAMlines = [samline.rstrip('\n') for samline in open(SAMfile)]
for thisline in SAMlines:
    if thisline[0]!="@":
        parts = thisline.split("\t")
        
        OligoID = parts[0]
        TrueMapping = parts[1]
        SoftClips = parts[5]
        Sequence = parts[9]
         
        if OligoID not in AllOligos.keys():
            NH = parts[11]
            NHparts = NH.split(":")
            AllOligos[OligoID] = NHparts[2]
            
            gc_perc = GCcontent(Sequence)
            GC_dict[OligoID] = gc_perc
            
            Matches[OligoID] = 0
            
        if OligoID not in Sequences.keys():
            if (TrueMapping=="0") | (TrueMapping=="256"):
                Sequences[OligoID] = Sequence
            
        for ThisMatch in re.findall('[0-9]*[A-Z]', SoftClips):
            if ThisMatch[-1:]=="M": # Add up all matching base pairs
                Matches[OligoID]=Matches[OligoID]+int(ThisMatch[:-1])
            elif (ThisMatch[-1:]=="D") | (ThisMatch[-1:]=="I"): # Penalty for insertions and deletions
                if OligoID not in Deletions.keys():
                    Deletions[OligoID]=int(ThisMatch[:-1])
                else:
                    Deletions[OligoID]=Deletions[OligoID]+int(ThisMatch[:-1])

# Calculate density score

DensityDict = {}
for ThisOligo in AllOligos.keys():
    TotalScore = 0
    if ThisOligo in Matches.keys():
        TotalScore = TotalScore+Matches[ThisOligo]
    if ThisOligo in Deletions.keys():
        TotalScore=TotalScore-Deletions[ThisOligo]
    Density = TotalScore/len(Sequences[ThisOligo]) # Normalise number of matches to length of oligo
    DensityDict[ThisOligo] = Density

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
TextOut.write("Chr\tStart\tStop\tFragment Start\tFragment Stop\tSide of fragment\tSequence\tTotal number of alignments\tDensity score\tRepeat length\tRepeat Class\tGC%\n")
for ThisOligo in AllOligos.keys():
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
        TextOut.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8:.2f}\t{9}\t{10}\t{11:.2f}\n".format(Chr,Start,Stop,FragStart,FragEnd,Side,Sequences[ThisOligo],AllOligos[ThisOligo],DensityDict[ThisOligo],RepeatLength,RepeatType,GC_dict[ThisOligo]))
        Written[OligoCoor] = 1
    #else:
        #print("Didn't write oligo {0}".format(OligoCoor)) ##Remove this later
TextOut.close()

sys.stdout.write("Finished running at: ")
TimeNow = datetime.datetime.now().time()
sys.stdout.write(str(TimeNow)+"\n")
