#!/usr/bin/env python

from Bio import SeqIO
import re,getopt,sys,os

# Functions
    
def usage():
    print("usage: FragExtract.py -b <bed file with coordinates> -g <genome build: hg18, hg19, mm9 or mm10> -e <restriction enzyme> -o <oligo size (bp)>")

input_file = ""
genome = ""
enzyme = ""
oligo_size = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'b:g:e:o:h',)
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
        elif opt == '-b':
            input_file = "../"+arg
        elif opt == '-g':
            genome = arg
        elif opt == '-e':
            enzyme = arg
        elif opt == '-o':
            oligo_size = arg
        else:
            usage()
            sys.exit(2)

Sequence_dict = {}
KillScript=0

if os.path.isfile(input_file)==False:
    print("Error: input file not recognised")
    KillScript=1

Organism = ""
if (genome=="hg18") | (genome=="hg19"):
    Organism="Homo_sapiens"
elif (genome=="mm9") | (genome=="mm10"):
    Organism="Mus_musculus"
else:
    print("Error: genome not recognised, please choose from hg18, hg19, mm9 or mm10")
    KillScript=1
    
Cut_sequence = ""
if enzyme=="DpnII":
    Cut_sequence="GATC"
elif enzyme=="NlaIII":
    Cut_sequence="CATG"
elif enzyme=="HindIII":
    Cut_sequence="AAGCTT"
else:
    print("Error: restriction enzyme is not recognised, please choose from DpnII, NlaIII or HindIII (case sensitive)")
    KillScript=1

try:
    oligo_value = int(oligo_size)
    if oligo_value<1:
        print("Oligo size invalid, it must be an integer greater than 0")
        KillScript=1
except ValueError:
    print("Error: oligo size invalid, it must be an integer greater than 0")
    KillScript=1
    
if KillScript==1:
    sys.exit(2)

Sequence = SeqIO.parse("/databank/igenomes/"+Organism+"/UCSC/"+genome+"/Sequence/WholeGenomeFasta/genome.fa", "fasta")
for seq_record in Sequence:
    Sequence_dict[seq_record.name] = seq_record.seq.upper()


FragmentCoordinates = {}    
def GenOligo(x,y,z):
    
    MoveAmount = len(Cut_sequence)-1
    
    SeqCheckUp = y-MoveAmount
    SeqCheckDown = z+MoveAmount
    Code = 0
    
    LeftUp=z-MoveAmount
    LeftDown=z
    LeftSeq = 0
    LeftUpCoor = 0
    LeftDownCoor = 0
    
    RightUp = y
    RightDown = y+MoveAmount
    RightSeq = 0
    RightUpCoor = 0
    RightDownCoor = 0
    
    FragmentLength = 0
    
    if re.search(Cut_sequence,str(Sequence_dict[Chr][SeqCheckUp:SeqCheckDown].upper())):
        Code = 1 # A code of 1 means check fragments either side
        LeftUpCoor=SeqCheckUp
        RightDownCoor=SeqCheckDown
    else:
        
        # Left oligo: walk left 1bp at a time until cut sequence or reach beginning of chromosome
        while bool(re.search(Cut_sequence,str(Sequence_dict[x][LeftUp:LeftDown])))==False:
            LeftUp=LeftUp-1
            if LeftUp<0:
                print("ran into start of chromosome")
                Code = 2 # A code of 2 means check fragment to the right
                break
            elif re.search(Cut_sequence,str(Sequence_dict[x][LeftUp:LeftDown])):
                LeftUpCoor = LeftUp
                LeftDownCoor = LeftUpCoor + oligo_value
                LeftSeq = Sequence_dict[x][LeftUpCoor:LeftDownCoor]
                Code = 4 # A code of 4 means fragment is fine
                break
 
        # Right oligo: walk right 1bp at a time until cut sequence or reach end of chromosome
        while bool(re.search(Cut_sequence,str(Sequence_dict[x][RightUp:RightDown])))==False:
            RightDown=RightDown+1
            if RightDown>len(Sequence_dict[x]):
                print("went over end of chromosome")
                Code = 3 # A code of 3 means check fragement to the left
                break
            elif re.search(Cut_sequence,str(Sequence_dict[x][RightUp:RightDown])):
                RightDownCoor = RightDown
                RightUpCoor = RightDownCoor-oligo_value
                RightSeq = Sequence_dict[x][RightUpCoor:RightDownCoor]
                if Code!=2:
                    Code = 4
                break
        
        # Check fragment length to exclude fragments < oligo size
        FragmentLength = RightDownCoor-LeftUpCoor
        if FragmentLength<oligo_value:
            if (Code!=2) & (Code!=3):
                Code = 1
        
    return Code,LeftSeq,LeftUpCoor,LeftDownCoor,RightSeq,RightUpCoor,RightDownCoor,FragmentLength

out_file = open("Report.txt","w")
error_file = open("ErrorReport.txt","w")
fasta_file = open("GeneratedOligos.fa","w")
bed_file = open("OligoCoor.bed","w")
gene_file = open("GeneAssociations.txt","w")
GroupCode = 1
gene_file.write("Fragment\tGenes\n")
Positions = [Position.rstrip('\n') for Position in open(input_file)]
for CurrentPosition in Positions:
    Chr,TSS1,TSS2,Gene = CurrentPosition.split('\t')
    
    PositionCoor = Chr+":"+TSS1+"-"+TSS2
    TSS1 = int(TSS1)
    TSS2 = int(TSS2)

    Code,LeftOligoSeq,LeftUpCoor,LeftDownCoor,RightOligoSeq,RightUpCoor,RightDownCoor,FragmentLength = GenOligo(Chr,TSS1,TSS2)
    if Code==4:
        FragCoor = Chr+":"+str(LeftUpCoor)+"-"+str(RightDownCoor)
        if FragCoor not in FragmentCoordinates.keys():
            gene_file.write(FragCoor+"\t"+Gene+"\n")
            out_file.write("Oligos were generated succesfully first time for {0}: L = {1}, R = {2}\n".format(PositionCoor,LeftOligoSeq,RightOligoSeq))
            FragmentCoordinates[FragCoor] = {}
            FragmentCoordinates[FragCoor].update({"Left": LeftOligoSeq})
            fasta_file.write(">"+Chr+":"+str(LeftUpCoor)+"-"+str(LeftDownCoor)+"-"+str(LeftUpCoor)+"-"+str(RightDownCoor)+"-L\n"+str(LeftOligoSeq)+"\n")
            FragmentCoordinates[FragCoor].update({"Right": RightOligoSeq})
            fasta_file.write(">"+Chr+":"+str(RightUpCoor)+"-"+str(RightDownCoor)+"-"+str(LeftUpCoor)+"-"+str(RightDownCoor)+"-R\n"+str(RightOligoSeq)+"\n")
            bed_file.write(Chr+"\t"+str(LeftUpCoor)+"\t"+str(LeftDownCoor)+"\n")
            bed_file.write(Chr+"\t"+str(RightUpCoor)+"\t"+str(RightDownCoor)+"\n")
        else:
            error_file.write("Position {0} was redundant with another position\n".format(PositionCoor))
            
    elif Code==1:
        error_file.write("Oligos could not be generated for {0} because the position was in cut site or fragment was too small\n".format(PositionCoor))
        LeftTSS1=LeftUpCoor-1
        LeftTSS2=LeftUpCoor
        RightTSS2=RightDownCoor+1
        RightTSS1=RightDownCoor
        
        Code_L,LeftOligoSeq_L,LeftUpCoor_L,LeftDownCoor_L,RightOligoSeq_L,RightUpCoor_L,RightDownCoor_L,FragmentLength_L = GenOligo(Chr,LeftTSS1,LeftTSS2)
        if Code_L==4:
            FragCoor = Chr+":"+str(LeftUpCoor_L)+"-"+str(RightDownCoor_L)
            if FragCoor not in FragmentCoordinates.keys():
                gene_file.write(FragCoor+"\t"+Gene+"\n")
                error_file.write("Oligos were generated succesfully on the fragment to the left\n")
                out_file.write("Oligos were generated succesfully on the adjancent left fragment for {0}: L = {1}, R = {2}\n".format(PositionCoor,LeftOligoSeq_L,RightOligoSeq_L))
                FragmentCoordinates[FragCoor] = {}
                FragmentCoordinates[FragCoor].update({"Left": LeftOligoSeq_L})
                fasta_file.write(">"+Chr+":"+str(LeftUpCoor_L)+"-"+str(LeftDownCoor_L)+"-"+str(LeftUpCoor_L)+"-"+str(RightDownCoor_L)+"-"+str(GroupCode)+"-L\n"+str(LeftOligoSeq_L)+"\n")
                FragmentCoordinates[FragCoor].update({"Right": RightOligoSeq_L})
                fasta_file.write(">"+Chr+":"+str(RightUpCoor_L)+"-"+str(RightDownCoor_L)+"-"+str(LeftUpCoor_L)+"-"+str(RightDownCoor_L)+"-"+str(GroupCode)+"-R\n"+str(RightOligoSeq_L)+"\n")
                bed_file.write(Chr+"\t"+str(LeftUpCoor_L)+"\t"+str(LeftDownCoor_L)+"\n")
                bed_file.write(Chr+"\t"+str(RightUpCoor_L)+"\t"+str(RightDownCoor_L)+"\n")
            else:
                error_file.write("Oligos were generated successfully for the adjacent left fragment but the fragment was redundant to another\n")  
        else:
            error_file.write("Oligos could not be generated in the left fragment\n")
        
        Code_R,LeftOligoSeq_R,LeftUpCoor_R,LeftDownCoor_R,RightOligoSeq_R,RightUpCoor_R,RightDownCoor_R,FragmentLength_R = GenOligo(Chr,RightTSS1,RightTSS2)
        if Code_R==4:
            FragCoor = Chr+":"+str(LeftUpCoor_R)+"-"+str(RightDownCoor_R)
            if FragCoor not in FragmentCoordinates.keys():
                gene_file.write(FragCoor+"\t"+Gene+"\n")
                error_file.write("Oligos were generated succesfully on the fragment to the right\n")
                out_file.write("Oligos were generated succesfully on the adjancent right fragment for {0}: L = {1}, R = {2}\n".format(PositionCoor,LeftOligoSeq_R,RightOligoSeq_R))
                FragmentCoordinates[FragCoor] = {}
                FragmentCoordinates[FragCoor].update({"Left": LeftOligoSeq_R})
                fasta_file.write(">"+Chr+":"+str(LeftUpCoor_R)+"-"+str(LeftDownCoor_R)+"-"+str(LeftUpCoor_R)+"-"+str(RightDownCoor_R)+"-"+str(GroupCode)+"-L\n"+str(LeftOligoSeq_R)+"\n")
                FragmentCoordinates[FragCoor].update({"Right": RightOligoSeq_R})
                fasta_file.write(">"+Chr+":"+str(RightUpCoor_R)+"-"+str(RightDownCoor_R)+"-"+str(LeftUpCoor_R)+"-"+str(RightDownCoor_R)+"-"+str(GroupCode)+"-R\n"+str(RightOligoSeq_R)+"\n")
                bed_file.write(Chr+"\t"+str(LeftUpCoor_R)+"\t"+str(LeftDownCoor_R)+"\n")
                bed_file.write(Chr+"\t"+str(RightUpCoor_R)+"\t"+str(RightDownCoor_R)+"\n")
            else:
                error_file.write("Oligos were generated successfully for the adjacent right fragment but the fragment was redundant to another\n")
        else:
            error_file.write("Oligos could not be generated in the right fragment\n")
        
        GroupCode=GroupCode+1
bed_file.close()
fasta_file.close()
out_file.close()
error_file.close()
gene_file.close()