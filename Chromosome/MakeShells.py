#!/usr/bin/env python

import getopt,sys

def usage():
    print("usage: MakeShells.py -g <genome build> -c <chromosome number> -u <upper limit>")
    
genome = ""
chromosome = ""
upperLimit = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], 'g:c:u:h',)
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
    out_file.write("/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ../chr"+chromosome+"_Oligos_"+Suffix+".fa --genomeDir /databank/igenomes/"+Organism+"/UCSC/"+genome+"/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr"+chromosome+"_"+Suffix+"_\n")
    #out_file.write("repeatmasker -noint -s -dir /t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/WholeChromosomes/mm9_chr11_DpnII_70bp/"+Suffix+"/ -species mouse ../chr11_Oligos_"+Suffix+".fa\n")
    out_file.write("repeatmasker -noint -s -dir ./ -species "+Species+" ../chr"+chromosome+"_Oligos_"+Suffix+".fa\n")
    #out_file.write("python /t1-home/nuffmed/jkerry/Python/STAR_depthGauge_v1-10.py -i "+Suffix+"\n")
    out_file.write("python /t1-data1/WTSA_Dev/jkerry/OligoDesign/Chromosome/OligoSTAR.py -i "+Suffix+" -c "+chromosome+"\n")
    #out_file.write("module load ucsctools\n")
    #out_file.write("wigToBigWig RM_plot_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt RM_plot_L_"+Suffix+".bw\nwigToBigWig RM_plot_R_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt RM_plot_R_"+Suffix+".bw\nwigToBigWig STAR_density_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_density_L_"+Suffix+".bw\nwigToBigWig STAR_density_R_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_density_R_"+Suffix+".bw\nwigToBigWig STAR_depthGauge_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_depthGauge_L_"+Suffix+".bw\nwigToBigWig STAR_depthGauge_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_depthGauge_L_"+Suffix+".bw\n")
    #out_file.write("mv *.bw /public/jkerry/\n")
    CounterLow = CounterLow+20000
    Counterhigh = Counterhigh+20000
out_file.close()

#CounterLow = 200001
#Counterhigh = 220000
#while Counterhigh<=740000:
#    Suffix = str(CounterLow)+"-"+str(Counterhigh)
#    out_file = open("density_"+Suffix+".sh","w")
#    out_file.write("mkdir "+Suffix+"\n")
#    out_file.write("cd "+Suffix+"\n")
#    out_file.write("/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_"+Suffix+".fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr1_"+Suffix+"_\n")
#    out_file.write("repeatmasker -noint -s -dir /t1-data1/WTSA_Dev/jkerry/CS_paper/"+Suffix+"/ -species mouse /t1-data1/WTSA_Dev/jkerry/CS_paper/chr1_Oligos_"+Suffix+".fa\n")
#    out_file.write("python /t1-home/nuffmed/jkerry/Python/STAR_density_v1-01.py -i "+Suffix+"\n")
#    CounterLow = CounterLow+20000
#    Counterhigh = Counterhigh+20000
#out_file.close()

#CounterLow = 1
#Counterhigh = 20000
#while Counterhigh<=740000:
#    Suffix = str(CounterLow)+"-"+str(Counterhigh)
#    out_file = open("density_"+Suffix+".sh","w")
#    out_file.write("cd "+Suffix+"\n")
#    out_file.write("python /t1-home/nuffmed/jkerry/Python/STAR_density_v1-10.py -i "+Suffix+"\n")
#    CounterLow = CounterLow+20000
#    Counterhigh = Counterhigh+20000
#out_file.close()