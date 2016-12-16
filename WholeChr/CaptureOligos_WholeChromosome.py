#!usr/bin/python

from Bio import SeqIO
import re,getopt,sys

def usage():
    print("usage: CaptureOligos_WholeChromosome.py -g <genome build> -c <chromosome number> -e <restriction enzyme> -o <oligo size (bp)> -r <region of chromsome (optional)>")
    print("For extended help: 'CaptureOligos_WholeChromosome.py -h'")
def help_info():
    print("usage: CaptureOligos_WholeChromosome.py -g <genome build> -c <chromosome number> -e <restriction enzyme> -o <oligo size (bp)>\n")
    print("-------------------------------------------------------------------------\n")
    print("CaptureOligos_WholeChromosome.py generates oligos adjacent to every restriction site of the supplied enzyme for an entire chromosome.\n")
    print("\tGenomes: choose from 'hg18', 'hg19', 'mm9' or 'mm10'\n")
    print("\tChromosome: supply the bare number or letter of the chromosome e.g. '7' or 'X'. Only one chromosome can be run at a time\n")
    print("\tRestriction enzymes: 'DpnII' (GATC), 'NlaIII' (CATG) or 'HindIII' (AAGCTT)\n")
    print("\tOligo size: specify the size of the oligos to be generated adjacent to restriction sites. Supply the number of base pairs e.g. '120'\n")
    print("Example (70bp oligos for DpnII on chr11 of mouse build mm9): 'CaptureOligos_WholeChromosome.py -g mm9 -c 11 -e DpnII -o 70'\n")
    
genome = ""
chromosome = ""
enzyme = ""
oligo_size = ""
region = ""

try:
    opts, args = getopt.getopt(sys.argv[1:], 'g:c:e:o:r:h',)
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
        elif opt == '-g':
            genome = arg
        elif opt == '-c':
            chromosome = arg
        elif opt == '-e':
            enzyme = arg
        elif opt == '-o':
            oligo_size = arg
        elif opt == '-r':
            region = arg
        else:
            usage()
            sys.exit(2)
            
## Check for errors

# Genome
            
KillScript = 0
MaxChromosome = 19
Organism = ""
if (genome=="hg18") | (genome=="hg19"):
    Organism="Homo_sapiens"
    MaxChromosome = 22
elif (genome=="mm9") | (genome=="mm10"):
    Organism="Mus_musculus"
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

# Restriction enzyme 
Cut_sequence = ""
if enzyme=="DpnII":
    Cut_sequence="GATC"
elif enzyme=="NlaIII":
    Cut_sequence="CATG"
elif enzyme=="HindIII":
    Cut_sequence="AAGCTT"
else:
    print("Restriction enzyme is not recognised, please choose from DpnII, NlaIII or HindIII (case sensitive)")
    KillScript=1
    
# Oligo size
    
try:
    oligo_value = int(oligo_size)
    if oligo_value<1:
        print("Oligo size invalid, it must be an integer greater than 0")
        KillScript=1
except ValueError:
    print("Oligo size invalid, it must be an integer greater than 0")
    KillScript=1
    
if KillScript==1:
    sys.exit(2)

## Store genome sequence in dictionary, by chromosome

Sequence_dict = {}
Sequence = SeqIO.parse("/databank/igenomes/"+Organism+"/UCSC/"+genome+"/Sequence/WholeGenomeFasta/genome.fa", "fasta")
for seq_record in Sequence:
    Sequence_dict[seq_record.name] = seq_record.seq.upper()

## Find positions of all restriction sites on the chromosome

Chr="chr"+chromosome
StartSeq = 1
StopSeq = len(Sequence_dict[Chr])
if region!="":    
    StartSeq,StopSeq = region.split("-")
StartSeq = int(StartSeq)-1    
StopSeq = int(StopSeq)
posList = []
p = re.compile(Cut_sequence)
for m in p.finditer(str(Sequence_dict[Chr][StartSeq:StopSeq])):
    posList.append(m.start()+StartSeq)

# Find all GATC sites

if region!="":
    fasta_file = open("Oligos_"+genome+"_chr"+chromosome+"_"+region+"_"+enzyme+"_"+oligo_size+"bp.fa","w")
else:
    fasta_file = open("Oligos_"+genome+"_chr"+chromosome+"_"+enzyme+"_"+oligo_size+"bp.fa","w")
ThisSite = 0
p = re.compile(Cut_sequence)
for m in p.finditer(str(Sequence_dict[Chr][StartSeq:StopSeq])):
    
    # Genereate oligo positions within sequence variable (70bp, including GATC site)
    #Position = m.start()
    Position = posList[ThisSite]
    
    LeftOligoStart = Position
    LeftOligoStop = Position + oligo_value
    if LeftOligoStop > StopSeq:
        LeftOligoStop=StopSeq
        
    ReadLeftStart = LeftOligoStart-StartSeq
    ReadLeftStop = LeftOligoStop-StartSeq
    
    LeftOligo = Sequence_dict[Chr][ReadLeftStart:ReadLeftStop]
    LeftCoor = Chr+":"+str(LeftOligoStart)+"-"+str(LeftOligoStop) # These coordinates match with the output from bedtools fastaFromBed but in WIG tracks they appear shifted 1bp to the left
    
    NextSite = ThisSite+1
    if NextSite<len(posList):
        NextPosition = posList[NextSite]
        RightOligoStart = NextPosition-(oligo_value-len(Cut_sequence))
        if RightOligoStart < StartSeq:
            RightOligoStart=StartSeq
        RightOligoStop = NextPosition+len(Cut_sequence)
        
        ReadRightStart = RightOligoStart-StartSeq
        ReadRightStop = RightOligoStop-StartSeq
    
        RightOligo = Sequence_dict[Chr][ReadRightStart:ReadRightStop]
        RightCoor = Chr+":"+str(RightOligoStart)+"-"+str(RightOligoStop) # These coordinates match with the output from bedtools fastaFromBed but in WIG tracks they appear shifted 1bp to the left
    
    # Check region is bigger than 70bp
    FragmentLeft=Position
    if NextSite<len(posList):
        FragmentLength=posList[NextSite]-posList[ThisSite]+len(Cut_sequence)
        FragmentRight=FragmentLeft+FragmentLength
        FragmentCoor=str(FragmentLeft)+"-"+str(FragmentRight)
        
        # fragment number = NextSite
        # fragment length = FragmentLength
        
        if FragmentLength>=oligo_value:
            if FragmentLength>=(oligo_value*2):
                fasta_file.write(">{0}-{1}-L\n{2}\n>{3}-{1}-R\n{4}\n".format(LeftCoor,FragmentCoor,LeftOligo,RightCoor,RightOligo))
            else:
                fasta_file.write(">{0}-{1}-L\n{2}\n".format(LeftCoor,FragmentCoor,LeftOligo))
    ThisSite=ThisSite+1
fasta_file.close()
