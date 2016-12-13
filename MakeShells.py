#!/usr/bin/env python
CounterLow = 1
Counterhigh = 20000
while Counterhigh<=460000:
    Suffix = str(CounterLow)+"-"+str(Counterhigh)
    out_file = open("Shell_"+Suffix+".sh","w")
    out_file.write("mkdir "+Suffix+"\n")
    out_file.write("cd "+Suffix+"\n")
    out_file.write("/package/rna-star/2.5.1b/bin/STAR --runThreadN 4 --readFilesIn ../chr11_Oligos_"+Suffix+".fa --genomeDir /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/STAR/ --genomeLoad NoSharedMemory --outFilterMultimapScoreRange 1000 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 110 --seedSearchStartLmax 4 --seedSearchLmax 20 --alignIntronMax 10 --seedPerWindowNmax 15 --seedMultimapNmax 11000 --winAnchorMultimapNmax 200 --limitOutSAMoneReadBytes 300000 --outFileNamePrefix chr11_"+Suffix+"_\n")
    out_file.write("repeatmasker -noint -s -dir /t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/WholeChromosomes/mm9_chr11_DpnII_70bp/"+Suffix+"/ -species mouse ../chr11_Oligos_"+Suffix+".fa\n")
    #out_file.write("python /t1-home/nuffmed/jkerry/Python/STAR_depthGauge_v1-10.py -i "+Suffix+"\n")
    out_file.write("python /t1-home/nuffmed/jkerry/Python/STAR_density_v1-10.py -i "+Suffix+"\n")
    out_file.write("module load ucsctools\n")
    out_file.write("wigToBigWig RM_plot_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt RM_plot_L_"+Suffix+".bw\nwigToBigWig RM_plot_R_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt RM_plot_R_"+Suffix+".bw\nwigToBigWig STAR_density_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_density_L_"+Suffix+".bw\nwigToBigWig STAR_density_R_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_density_R_"+Suffix+".bw\nwigToBigWig STAR_depthGauge_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_depthGauge_L_"+Suffix+".bw\nwigToBigWig STAR_depthGauge_L_"+Suffix+".wig /databank/igenomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/chr_sizes.txt STAR_depthGauge_L_"+Suffix+".bw\n")
    out_file.write("mv *.bw /public/jkerry/\n")
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