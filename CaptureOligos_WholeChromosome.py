#!usr/bin/python

from Bio import SeqIO
import re,getopt,sys

def usage():
    print("usage: CaptureOligos_WholeChromosome.py -g <genome build> -c <chromosome number> -e <restriction enzyme> -o <oligo size (bp)>")
    print("for extended help, type CaptureOligos_WholeChromosome.py -h")
def help():
    print("usage: CaptureOligos_WholeChromosome.py -g <genome build> -c <chromosome number> -e <restriction enzyme> -o <oligo size (bp)>")
    print("CaptureOligos_WholeChromosome.py generates oligos adjacent to every restriction site of the supplied enzyme for an entire chromosome.")
    print("Genomes: choose from 'hg18', 'hg19', 'mm9' or 'mm10'")
    print("Chromosome: supply the bare number of letter e.g. '7' or 'X'. Only one chromosome can be run at a time")
    print("Oligo size: specify the size of the oligos to be generated adjacent to restriction sites. Supply the number of base pairs e.g. '70'")
    
genome = ""
chromosome = ""
enzyme = ""
oligo_size = ""

try:
    opts, args = getopt.getopt(sys.argv[1:], 'g:c:e:o:h',)
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
        elif opt == '-e':
            enzyme = arg
        elif opt == '-o':
            oligo_size = arg
        else:
            usage()
            sys.exit(2)


# Put sequence of chr1 into a single variable
#chr1_file = "/t1-data1/WTSA_Dev/jkerry/CS_paper/Mus_musculus.NCBIM37.67.dna.chromosome.1.fa"
chr1_file = "/t1-data1/WTSA_Dev/jkerry/CS_paper/chr1.fa"
fasta_file = open("chr1_DpnII_70bpOligos.fa","w")
Complete_seq = ""
chr1_lines = [chr1_line.rstrip('\n') for chr1_line in open(chr1_file)]
for thisLine in chr1_lines[1:]:
    Complete_seq = Complete_seq+""+thisLine.upper()

# re test
#p = re.compile("GATC")
#testSeq = "TAGCGAGAGTTGCAGGATCAGATGCGAGGTGCGAGTGCGATGGATGCGATCCAAAGCTGAGA"
#for m in p.finditer(testSeq):
#    print m.start(), m.group()
#    Position = m.start()
#    StartOligo = Position-10
#    LeftOligo = testSeq[StartOligo:Position]
#    print LeftOligo

# Log positions of GATC in array
Chr="chr1"
StartSeq = 1
StartSeq = StartSeq-1
StopSeq = 197195432
posList = []
p = re.compile("GATC")
for m in p.finditer(Complete_seq):
    posList.append(m.start()+StartSeq)
#print posList

# Find all GATC sites

ThisSite = 0
p = re.compile("GATC")
for m in p.finditer(Complete_seq):
    
    # Genereate oligo positions within sequence variable (70bp, including GATC site)
    #Position = m.start()
    Position = posList[ThisSite]
    
    LeftOligoStart = Position
    LeftOligoStop = Position + 70
    if LeftOligoStop > StopSeq:
        LeftOligoStop=StopSeq
        
    ReadLeftStart = LeftOligoStart-StartSeq
    ReadLeftStop = LeftOligoStop-StartSeq
    
    LeftOligo = Complete_seq[ReadLeftStart:ReadLeftStop]
    LeftCoor = Chr+":"+str(LeftOligoStart)+"-"+str(LeftOligoStop) # These coordinates match with the output from bedtools fastaFromBed but in WIG tracks they appear shifted 1bp to the left
    
    NextSite = ThisSite+1
    if NextSite<len(posList):
        NextPosition = posList[NextSite]
        RightOligoStart = NextPosition-66
        if RightOligoStart < StartSeq:
            RightOligoStart=StartSeq
        RightOligoStop = NextPosition+4
        
        ReadRightStart = RightOligoStart-StartSeq
        ReadRightStop = RightOligoStop-StartSeq
    
        RightOligo = Complete_seq[ReadRightStart:ReadRightStop]
        RightCoor = Chr+":"+str(RightOligoStart)+"-"+str(RightOligoStop) # These coordinates match with the output from bedtools fastaFromBed but in WIG tracks they appear shifted 1bp to the left
    
    # Check region is bigger than 70bp
    FragmentLeft=Position
    if NextSite<len(posList):
        FragmentLength=posList[NextSite]-posList[ThisSite]+4
        FragmentRight=FragmentLeft+FragmentLength
        FragmentCoor=str(FragmentLeft)+"-"+str(FragmentRight)
        
        # fragment number = NextSite
        # fragment length = FragmentLength
        
        if FragmentLength>=70:
            if FragmentLength>=140:
                fasta_file.write(">{0}-{1}-L\n{2}\n>{3}-{1}-R\n{4}\n".format(LeftCoor,FragmentCoor,LeftOligo,RightCoor,RightOligo))
            else:
                fasta_file.write(">{0}-{1}-L\n{2}\n".format(LeftCoor,FragmentCoor,LeftOligo))
    ThisSite=ThisSite+1
fasta_file.close()
