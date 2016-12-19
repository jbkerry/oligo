#!usr/bin/perl

use Math::Round;


my %depthgauge;
my %Maffed;
my %used;
my %oligoCoverage;
my %OligoValue;
my %OligoSize;
my %StoredSeq;
my %dup_counts;
my %ProblemSequences;
my %qGapBases;

my $MAfLine = 0;
my $MafHeader;
my $MafSeq;


#my $MafFile = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/GeneratedOligos.fa";
my $MafFile = "./GeneratedOligos.fa";

open (MAF, $MafFile);	
while (<MAF>)
{
    chomp;	    
    $MAfLine++;
	
    if ($MAfLine == 1)
    {
	$MafHeader = $_;
	$MafHeader =~ s/>//g;

	my ($Mchr, $Mstart, $Mstop, $Mfragstart, $Mfragstop, $Mside) = split(/[\:\-]/, $MafHeader);
	my $Mcounter = $Mstart;
	$Maffed{$MafHeader}=0;
	until ($Mcounter == $Mstop + 1)
	{
	    $depthgauge{$Mchr}{$Mcounter}=0;
	    $oligoCoverage{$Mchr}{$Mcounter}++;
	    $Mcounter++;
	}
    }
	
    if ($MAfLine == 2)
    {
	$MafSeq = $_;
	$StoredSeq{$MafHeader}=$MafSeq;
	$MAfLine=0;
    }
}	
close MAF;


#my $file = "/t1-data1/WTSA_Dev/jkerry/CaptureC/WholeGenome/Promoters/RefSeq_data/Dec2016/GeneratedOligos.pslx";
my $file = "./GeneratedOligos.pslx";
my $linenumber = 0;
open (BLAT, $file);	
while (<BLAT>)
{
    chomp;
    $linenumber++;
    if ($linenumber <= 5) { next; }
	
    my ($match, $mismatch, $repmatch, $ns, $qgap, $qgapbases, $tgap, $tgapbases, $strand, $query, $qsize, $qstart, $qend, $tname, $tsize, $tstart, $tend, $blockcount, $blocksize, $qtarts, $tstarts, $matchSeq)=split(/\s+/);
    if ($tname =~ /_/g){next;}
	
    my $percent = ($qsize / 100);	
    my $percent_match = ($match / $percent);	
    if ($percent_match >= 70) { $dup_counts{$query}++; }
    
    $qGapBases{$query}+=$qgapbases;
    my ($chr, $start, $stop, $fragstart, $fragstop, $side)=split(/[\:\-]/, $query);		
    #if ($linenumber<=100) {
    #    print $chr."\t".$start."\t".$stop."\t".$rest."\n";
    #}
        
    unless (exists $used{$query})
    {
	my $counter = $start;			
	until ($counter == ($start + $qsize + 1))
	{
	    $oligoCoverage{$chr}{$counter}++;
	    $counter++; 
	} 
    }
	
    $Maffed{$query}++;
    $used{$query}++;
    $OligoValue{$query} = 0;		
	
    ### This is where it tots up all blat hits for a given oligo
	
    my $loopstart = $start + $qstart;
    my $counter = $qstart;
	
    until ($counter == ($qend + 1))
    {
	$depthgauge{$chr}{$loopstart}++;
	$counter++;
	$loopstart++;		
    }
}
close BLAT;

foreach my $StoredID (keys %Maffed)
{	
    my ($storedchr, $storedstr, $storedstp, $storedFragstr, $storedFragstp, $storedSide)=split(/[\:\-]/, $StoredID); 
    my $posCounter = $storedstr;

    my $size = $storedstp - $storedstr;
    $OligoSize{$StoredID}=$size;

    until ($posCounter == $storedstp + 1)
    {
	if (exists $depthgauge{$storedchr}{$posCounter})
	{
	    $OligoValue{$StoredID}+=$depthgauge{$storedchr}{$posCounter};
	    $posCounter++;
	}
	else
	{
	    print "Warning: entry missing from blat file\n";
	}
    }
    $OligoValue{$StoredID}-=$qGapBases{$StoredID};
}

my %Groups;
foreach my $Did (sort keys %OligoValue)
{
    my $density = $OligoValue{$Did}/$OligoSize{$Did};
    my $density_round = nearest(.01, $density);
    my @testArray=split(/[\:\-]/, $Did);
    my $length = scalar(@testArray);
    if ($length==7) {
	my ($Dchr, $Dstr, $Dstp, $Dfragstr, $Dfragstp, $Dgroup, $Dside)=split(/[\:\-]/, $Did);
	my $FragCoor = $Dchr.":".$Dfragstr."-".$Dfragstp;
	#print "density= ".$Dgroup."\n";
	if (exists $Groups{$Dgroup}{$FragCoor}) {
	    $Groups{$Dgroup}{$FragCoor}+=$density_round;
	}
	else {
	   $Groups{$Dgroup}{$FragCoor}=$density_round;
	}
	#if (not exists $Groups{$Dgroup}) {
	#    $Groups{$Dgroup}=$density_round;
	#}
	#elsif ($density_round<$Groups{$Dgroup}) {
	#    $Groups{$Dgroup}=$density_round;
	#}
    }
}

my %LowestGroup;
foreach my $key (sort {$a <=> $b} keys %Groups) {
    #print "key = ".$key."\n";
    foreach my $subkey (keys $Groups{$key}) {
	#print "key =".$key.", subkey = ".$subkey.", density = ".$Groups{$key}{$subkey}."\n";
	if (exists $LowestGroup{$key}) {
	    if ($Groups{$key}{$subkey}<$LowestGroup{$key}) {
		$LowestGroup{$key}=$Groups{$key}{$subkey};
	    } 
	}
	else {
	    $LowestGroup{$key}=$Groups{$key}{$subkey}
	}
    }
}


## Repeat masker

my %LengthHash;
my %RMscoreHash = ();
my %SSRlengthHash = ();
open(RMout, "GeneratedOligos.fa.out") or die $!;
<RMout>;
<RMout>;
<RMout>;
while (<RMout>) {
    my $RepeatLine = $_;
    chomp($RepeatLine);
    (my $whiteSpace, my $sw_score, my $perc_div, my $perc_del, my $perc_ins, my $Qname, my $Qstart, my $Qstop, my $Qleft, my $PlusSign, my $Rname, my $Rclass, my $Rstart, my $Rstop, my $Rleft, my $LineID) = split(/\s+/, $RepeatLine);
    
    (my $OligoChr, my $OligoStart, my $OligoStop, my $FragStart, my $FragStop, my $Group, my $Side) = "";
    my @testArray=split(/[\:\-]/, $Qname);
    my $length = scalar(@testArray);
    if ($length==7) {
	($OligoChr, $OligoStart, $OligoStop, $FragStart, $FragStop, $Group, $Side) = split(/[\:\-]/, $Qname);
    }
    else {
	($OligoChr, $OligoStart, $OligoStop, $FragStart, $FragStop, $Side) = split(/[\:\-]/, $Qname);
    }
    
    my $SSRlength = $Qstop - $Qstart;
    $SSRlengthHash{$Qname} = $SSRlength;
    
    my $coor = $OligoChr.":".$OligoStart."-".$OligoStop;
    if (exists $LengthHash{$coor}) {
        if ($SSRlength>$LengthHash{$coor}) {
            $LengthHash{$coor}=$SSRlength
        }
    }
    else {
        $LengthHash{$coor}=$SSRlength;
    }
    #print "coor = ".$coor.", repeat = ".$SSRlength."\n";
    #if ($SSRlength>=15 && $SSRlength<=30 && $DensityHash{$coor}<=30) {
    #    print $SSRlength."\n";
    #}
}
close(RMout);

my %Written;
open(OUTPUT, ">Oligos_filtered.txt") or die $!;
print OUTPUT "Chr\tOligo Start\tOligo Stop\tFragment Start\tFragment Stop\tSide of Fragment\tDensity\tRepeat Length\tSequence\n";
foreach my $Did (sort keys %OligoValue)
{
    
    my $density = $OligoValue{$Did}/$OligoSize{$Did};
    my $density_round = nearest(.01, $density);
    my $Write=0;
    my ($Dchr, $Dstr, $Dstp, $Dfragstr, $Dfragstp, $Dgroup, $Dside) = "";
    my $OligoCoor = "";
    my @testArray=split(/[\:\-]/, $Did);
    my $length = scalar(@testArray);
    if ($length==7) {
	($Dchr, $Dstr, $Dstp, $Dfragstr, $Dfragstp, $Dgroup, $Dside)=split(/[\:\-]/, $Did);
	#print "This equals, ".$Dgroup."\t".$Dside."\n";
	$OligoCoor = $Dchr.":".$Dstr."-".$Dstp;
	my $FragCoor = $Dchr.":".$Dfragstr."-".$Dfragstp;
	unless (exists $Written{$OligoCoor}) {
	    $Write=1;
	    if (exists $LowestGroup{$Dgroup}) {
		if ($Groups{$Dgroup}{$FragCoor}==$LowestGroup{$Dgroup}) {
		    $Write=1;
		}
		else {
		    $Write=0;
		}
	    } 
	}
    }
    else {
	($Dchr, $Dstr, $Dstp, $Dfragstr, $Dfragstp, $Dside)=split(/[\:\-]/, $Did);
	$OligoCoor = $Dchr.":".$Dstr."-".$Dstp;
	unless (exists $Written{$OligoCoor}) {
	    $Write=1;
	}
    }
	
    if ($Write==1) {
	unless ($density_round>30 || $LengthHash{$OligoCoor}>30) {
	    my $PrintRepeats = "NA";
	    if ($LengthHash{$OligoCoor}!="") {
		$PrintRepeats = $LengthHash{$OligoCoor};
	    }
	    #print OUTPUT $Dchr.":".$Dstr."-".$Dstp."-".$Dfragstr."-".$Dfragstp."-".$Dside."\t".$density_round."\t".$PrintRepeats."\t".$StoredSeq{$Did}."\n";
	    print OUTPUT $Dchr."\t".$Dstr."\t".$Dstp."\t".$Dfragstr."\t".$Dfragstp."\t".$Dside."\t".$density_round."\t".$PrintRepeats."\t".$StoredSeq{$Did}."\n";
	    $Written{$OligoCoor}=1;
	}
    }
}
close(OUTPUT);
