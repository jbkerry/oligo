#!usr/bin/python

from __future__ import division
import re

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

MafFile = "./GeneratedOligos.fa"

MafLines = [MafLineRead.rstrip('\n') for MafLineRead in open(MafFile)]
for ThisMafLine in MafLines:    
    MafLine=MafLine+1
    if MafLine==1:
        MafHeader = ThisMafLine
        MafHeader = re.sub('>','',MafHeader)

        Mchr, Mstart, Mstop, Mfragstart, Mfragstop, Mside = re.split("\W+",MafHeader,5)
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
        MafLine=0       


BLATfile = "./GeneratedOligos.psl";
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
    Oligochr, start, stop, fragstart, fragstop, side=re.split("\W+",query,5)            
        
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
    storedchr, storedstr, storedstp, storedFragstr, storedFragstp, storedSide=re.split("\W+",StoredID,5); 
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

Groups = {}
for Did in sorted(OligoValue.keys()):
    density = OligoValue[Did]/OligoSize[Did];
    density_round = float("{0:.2f}".format(density))
    if len(re.split("\W+",Did))==7:
        Dchr, Dstr, Dstp, Dfragstr, Dfragstp, Dgroup, Dside=re.split("\W+",Did)
        if Dgroup not in Groups.keys():
            Groups[Dgroup] = {}
        FragCoor = Dchr+":"+Dfragstr+"-"+Dfragstp
        if FragCoor in Groups[Dgroup].keys():
            Groups[Dgroup][FragCoor]=Groups[Dgroup][FragCoor]+density_round
        else:
            Groups[Dgroup].update({FragCoor: density_round})
       #if (not exists $Groups{$Dgroup}) {
       #    $Groups{$Dgroup}=$density_round;
       #}
       #elsif ($density_round<$Groups{$Dgroup}) {
       #    $Groups{$Dgroup}=$density_round;
       #}

LowestGroup = {}
for key in sorted(Groups.keys()):
    for subkey in Groups[key].keys():
        if key in LowestGroup.keys():
            if Groups[key][subkey]<LowestGroup[key]:
                LowestGroup[key]=Groups[key][subkey]
        else:
            LowestGroup[key]=Groups[key][subkey]

## Repeat masker

LengthDict = {}
RMscoreDict = {}
SSRlengthDict = {}
RM_file="GeneratedOligos.fa.out"
RMLines = [RMLine.rstrip('\n') for RMLine in open(RM_file)]
for ThisRMLine in RMLines[3:]:
    RepeatLine = ThisRMLine
    whiteSpace, sw_score, perc_div, perc_del, perc_ins, Qname, Qstart, Qstop, Qleft, PlusSign, Rname, Rclass, Rstart, Rstop, Rleft, LineID, whiteSpace2 = re.split("\s+",RepeatLine)
    
    #OligoChr = OligoStart = OligoStop = FragStart = FragStop = Group = Side = "";
    if len(re.split("\W+",Qname))==7:
       OligoChr, OligoStart, OligoStop, FragStart, FragStop, Group, Side = re.split("\W+", Qname)
    else:
       OligoChr, OligoStart, OligoStop, FragStart, FragStop, Side = re.split("\W+", Qname)
    
    SSRlength = int(Qstop) - int(Qstart)
    SSRlengthDict[Qname] = SSRlength
    
    coor = OligoChr+":"+OligoStart+"-"+OligoStop
    if coor in LengthDict.keys():
        if (SSRlength>LengthDict[coor]):
            LengthDict[coor]=SSRlength
    else:
        LengthDict[coor]=SSRlength
    #print "coor = ".$coor.", repeat = ".$SSRlength."\n";
    #if ($SSRlength>=15 && $SSRlength<=30 && $DensityHash{$coor}<=30) {
    #    print $SSRlength."\n";
    #}
    
    

Written = {}
OUTPUT = open("Oligos_filtered.txt","w")
OUTPUT.write("Chr\tOligo Start\tOligo Stop\tFragment Start\tFragment Stop\tSide of Fragment\tDensity\tRepeat Length\tSequence\n")
for Did in sorted(OligoValue.keys()):
    density = OligoValue[Did]/OligoSize[Did]
    density_round = float("{0:.2f}".format(density))
    Write=0
    #$Dchr, $Dstr, $Dstp, $Dfragstr, $Dfragstp, $Dgroup, $Dside) = "";
    #OligoCoor = "";
    if len(re.split("\W+",Did))==7:
        Dchr, Dstr, Dstp, Dfragstr, Dfragstp, Dgroup, Dside=re.split("\W+", Did)
        #print "This equals, ".$Dgroup."\t".$Dside."\n";
        OligoCoor = Dchr+":"+Dstr+"-"+Dstp
        FragCoor = Dchr+":"+Dfragstr+"-"+Dfragstp
        if OligoCoor not in Written.keys():
            Write=1
            if  Dgroup in LowestGroup.keys():
                if Groups[Dgroup][FragCoor]==LowestGroup[Dgroup]:
                    Write=1
                else:
                    Write=0
    else:
        Dchr, Dstr, Dstp, Dfragstr, Dfragstp, Dside=re.split("\W+", Did)
        OligoCoor = Dchr+":"+Dstr+"-"+Dstp
        if OligoCoor not in Written.keys():
            Write=1
       
    if Write==1:
        if OligoCoor not in LengthDict.keys():
            LengthDict[OligoCoor]=0
        if (density_round<=30) & (LengthDict[OligoCoor]<=30):
            #print OUTPUT $Dchr.":".$Dstr."-".$Dstp."-".$Dfragstr."-".$Dfragstp."-".$Dside."\t".$density_round."\t".$PrintRepeats."\t".$StoredSeq{$Did}."\n";
            OUTPUT.write(Dchr+"\t"+Dstr+"\t"+Dstp+"\t"+Dfragstr+"\t"+Dfragstp+"\t"+Dside+"\t"+str(density_round)+"\t"+str(LengthDict[OligoCoor])+"\t"+StoredSeq[Did]+"\n")
            Written[OligoCoor]=1
