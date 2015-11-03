#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($help, $tabfile, $fastafile, $fastaRoot, $emblOut, $gffOut, $alnOut, @cmsearchOptions, $modelsDir, $model, $modelRoot, $cmsearchThresh, $cmsearchEvalue, $genome, $gene, $verbose, @warnings);
my $cmsearchThreshDefaultGenome = 16; #bits score threshold
my $modelGenome                 = 'genome.cm';
my $cmsearchThreshDefaultGene   = 14; #bits score threshold
my $modelGene                   = 'gene.cm';

&GetOptions( 
    "f|fastafile=s"     => \$fastafile,
    "p|prefix=s"        => \$fastaRoot,
    "t|tabfile=s"       => \$tabfile,
    "md|modeldir=s"     => \$modelsDir,
    "m|model=s"         => \$model,
    "th|thresh=s"       => \$cmsearchThresh,
    "e|evalue=s"        => \$cmsearchEvalue,
    "o|infOption=s@"    => \@cmsearchOptions,
    "specific|genome"   => \$genome,
    "sensitive|gene"    => \$gene,
    "e|embl"            => \$emblOut,
    "g|gff"             => \$gffOut,
    "a|aln"             => \$alnOut,
    "v|verbose"         => \$verbose,
    "h|help"            => \$help 
    );

if( $help ) {
    &help();
    exit(1);
}
elsif ((not defined $fastafile) && (not defined $tabfile)){
    print "ERROR: need a fastafile or a tabfile to work with!\n";
    &help();
    exit(1);    
}
elsif ((defined $fastafile) && ((not defined $modelsDir) && (not defined $ENV{'RNIE'}) && (not defined $model))){
    print "ERROR: neither model dir nor RNIE environment variable nor a model has been set!\n";
    &help();
    exit(1);    
}
elsif ((defined $genome) && (defined $gene)){
    print "ERROR: Idiot! Both genome and gene options have been set!\n";
    &help();
    exit(1);    
}

$modelsDir  = $ENV{'RNIE'} . '/models'           if ((not defined $modelsDir)  && (defined $ENV{'RNIE'}) && (length($ENV{'RNIE'})>0));
if ((not defined $gene) or (defined $genome)){
    $genome         = 1;
    $model          = $modelGenome                 if ( not defined $model);    
    $cmsearchThresh = $cmsearchThreshDefaultGenome if ( not defined $cmsearchThresh && not defined $cmsearchEvalue);
}
else {
    $gene           = 1;
    $model          = $modelGene                 if ( not defined $model);
    $cmsearchThresh = $cmsearchThreshDefaultGene if ( not defined $cmsearchThresh && not defined $cmsearchEvalue); 
}
$fastaRoot  = fileRoot($fastafile)               if ((not defined $fastaRoot) && (defined $fastafile) && (length($fastafile)>0));

$gffOut     = 1                                  if ((not defined $gffOut)     && (not defined $emblOut) && (not defined $alnOut));
$modelRoot  = fileRoot($model)                   if (     defined $model);
$fastaRoot .= '-' . $modelRoot                   if (     defined $modelRoot && defined $fastaRoot);
$fastaRoot  =~ s/\.cm$//                         if (     defined $fastaRoot);      
$fastaRoot .= 'Mode'                             if (     defined $fastaRoot);

printf "Config:\n"                               if $verbose;
print  "\tRNIE directory:  [".$ENV{'RNIE'}."]\n" if ($verbose && defined $ENV{'RNIE'});
printf "\tfastaFile input: [$fastafile]\n"       if ($verbose && defined $fastafile);
printf "\ttabFile input:   [$tabfile]\n"         if ($verbose && defined $tabfile);
printf "\tfastaRoot:       [$fastaRoot]\n"       if ($verbose && defined $fastaRoot);
printf "\tmodelsDir:       [$modelsDir]\n"       if ($verbose && defined $modelsDir);
printf "\tgenome:          [$genome]\n"          if ($verbose && defined $genome);
printf "\tgene:            [$gene]\n"            if ($verbose && defined $gene);
printf "\tmodel:           [$model]\n"           if ($verbose && defined $model);
printf "\tgffOut:          [$gffOut]\n"          if ($verbose && defined $gffOut);
printf "\temblOut:         [$emblOut]\n"         if ($verbose && defined $emblOut);
printf "\talnOut:          [$alnOut]\n"          if ($verbose && defined $alnOut);
printf "\tcmsearchThresh:  [$cmsearchThresh]\n"  if ($verbose && defined $cmsearchThresh);
printf "\tcmsearchEvalue:  [$cmsearchEvalue]\n"  if ($verbose && defined $cmsearchEvalue);

# make sure files are writable by group
umask(002);

if(!$tabfile && $fastafile){
    $tabfile = $fastaRoot . '.tabfile';
    my $cmsearchOut = $fastaRoot . '.cmsearch';
    my $cmsearchOptions;
    my $cmsearchThresholdString;
    if(defined $cmsearchEvalue){
	$cmsearchThresholdString = "-E $cmsearchEvalue";
    }
    else {
	$cmsearchThresholdString = "-T $cmsearchThresh";
    }
    
    $cmsearchOptions = "  $cmsearchThresholdString -g --fil-no-qdb --fil-T-hmm 2 --cyk --beta 0.05 " if defined $genome;
    $cmsearchOptions = "  $cmsearchThresholdString -g --fil-no-qdb --fil-no-hmm --no-qdb --inside  " if defined $gene;
    my $modelString = '';
    $modelsDir =~ s/\/$// if (defined $modelsDir);
    if (defined $model && -s $model){
	$modelString = $model;
    }
    elsif (defined $modelsDir && defined $model && (-s $modelsDir . '/' . $model)){
	$modelString = $modelsDir . '/' . $model;
    }
    elsif (defined $ENV{'RNIE'} && defined $model && (-s $ENV{'RNIE'} . '/' . $model)){
	$modelString = $ENV{'RNIE'} . '/' . $model;
    }
    else {
	print "ERROR: failed to find a model file!\n";
	print "       I check [$model], [" . $modelsDir . '/' . $model . "] & [" . $ENV{'RNIE'} . '/' . $model . "]\n";
	print "       you can give a full path to your model file with -m <str>, set a model directory with -md <dir> and a model file with -m <filename>\n";
	print "       or set the RNIE environment variable use one of the default models.\n";
	&help();
	exit(1);        
    }
    

    $cmsearchOptions .= join(' ', @cmsearchOptions) if @cmsearchOptions;
    my $cmsearchCmd = "cmsearch $cmsearchOptions --tabfile $tabfile  $modelString $fastafile > $cmsearchOut";
    print "Running:\n[$cmsearchCmd]\n" if $verbose;
    system($cmsearchCmd) and die "FATAL: failed to execute [$cmsearchCmd].\n[$!]";
}

my @cmNames;
if(-s $tabfile){
    #print "tabfile: [$tabfile]\n" if (defined $verbose);
    $fastaRoot = fileRoot($tabfile) if ((not defined $fastaRoot) && (defined $tabfile));
    #print "tabfileRoot: [$fastaRoot]\n" if (defined $verbose);
    my (%store);
    my ($counter,$cmName)=(10000,'');
#    open(F, "grep -v ^'#' $tabfile | sort -k6nr | ") or die "FATAL: could not open pipe for reading $tabfile\n[$!]";
    open(F, "< $tabfile") or die "FATAL: could not open $tabfile for reading\n[$!]";
  TABFILE: while (my $tabline = <F>){
      my( $model, $name, $start, $end, $qStart,$qEnd, $bits, $evalue, $gc );
      if($tabline=~/^# CM:\s+(\S+)/){
	  $cmName=$1;
	  push(@cmNames, $cmName);
      }
      elsif( (( $model, $name, $start, $end, $qStart,$qEnd, $bits, $evalue, $gc ) =
	    ($tabline=~/^\s*(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/)) || 
	            (( $name, $start, $end, $qStart,$qEnd, $bits, $evalue, $gc ) =
	            ($tabline=~/^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)$/)) ) {#Deals with both 1.0.0 and 1.0.2 formats!


	  next TABFILE if (defined $cmsearchThresh && $bits < $cmsearchThresh);
	  my $strand = checkStrand($start,$end);
	  ###################################
	  #Overlap checks. We only want the highest scoring instance
	  #of each overlapping frag.
	  my $overlaps = overlapCheck($name,$start,$end,$strand,$bits, $evalue,$cmName,\%store);
	  next TABFILE if $overlaps; 
	  
	  my $tag = 'terminator' . $counter;
	  push( @{ $store{$name} }, { 'start'   => $start,
				      'end'     => $end,
				      'strand'  => $strand,
				      'bits'    => $bits,
				      'evalue'  => $evalue,
				      'sumBits' => $bits,
				      'tag'     => $tag,
				      'cmName'  => $cmName,
				      'cmNames' => $cmName . "=" . $bits
		} );
	  
	  $counter++;
      }
  }
    
    #Print results to output file(s)
    print "\tGFF output:      [$fastaRoot\-rnie.gff]\n"  if ($verbose && $gffOut);
    open(GFF, "> $fastaRoot\-rnie.gff")  if $gffOut;
    print "\tEMBL output:     [$fastaRoot\-rnie.embl]\n" if ($verbose && $emblOut);
    open(EMBL, "> $fastaRoot\-rnie.embl") if $emblOut;
    my (%alnFileNames,%alnFilePointers);
    if ($alnOut){
	foreach my $cmName (@cmNames){
	    my $alnFileName = "$fastaRoot\-rnie.$cmName\.namefile";
	    $alnFileNames{$cmName} = $alnFileName;
	    open($alnFilePointers{$cmName}, "> $alnFileName");
	}
	
	if (not -s $fastafile . '.ssi'){
	    system("esl-sfetch --index $fastafile ") 
		and die "FATAL: failed to execute esl-sfetch.\n[$!]";
	}
    }

    foreach my $name ( sort{ $a cmp $b } keys %store) {
	my $nameSuffix=';';
	$nameSuffix='' if $name !~ /\;$/;
	printf EMBL "ID   $name$nameSuffix\n" if $emblOut;
	foreach my $feature (@{ $store{$name} }){
	    my ($start, $end, $bits, $evalue, $sumBits, $tag, $strand)=($feature->{'start'}, $feature->{'end'}, $feature->{'bits'}, $feature->{'evalue'}, $feature->{'sumBits'}, $feature->{'tag'}, $feature->{'strand'});
	    
	    if ($emblOut){
		my $coords=$start . ".." . $end; 
		$coords="complement(" . $start . ".." . $end . ")" if $start > $end;
		printf EMBL "FT   terminator       $coords
FT                    /note=\42Predicted Rho independent terminator using cmsearch\42
FT                    /note=\42Score " . $bits . "\42
FT                    /note=\42E value " . $evalue . "\42
FT                    /gene=\42terminator" . $tag . "\42\n";
	    }
	    
	    if ($gffOut){
		my $strandString = '+';
		$strandString = '-' if $strand == -1;
		my ($gffStart, $gffEnd)=($start,$end);
		($gffStart, $gffEnd)=($end,$start) if $gffStart>$gffEnd;
		printf GFF "$name\tRNIE\tterminator\t$gffStart\t$gffEnd\t$bits\t$strandString\t\.\tID=$tag;E_value=$evalue;Bits=$bits;sumBits=$sumBits;Note=Predicted Rho independent terminator using cmsearch,%s\n", $feature->{'cmNames'};
	    }
	    
	    if ($alnOut){
		printf {$alnFilePointers{$feature->{'cmName'}}} "$name/$start\-$end\t$start\t$end\t$name\n";
	    }
	}
    }
    close(EMBL) if $emblOut;
    close(GFF)  if $gffOut;
    if ($alnOut){
	foreach my $cmName (@cmNames){
	    close($alnFilePointers{$cmName});
	    system("esl-sfetch -o ". $alnFileNames{$cmName} .".fa -C -f $fastafile " . $alnFileNames{$cmName}) 
		and die "FATAL: failed to execute esl-sfetch.\n[$!]";
	    my $outFile = "$fastaRoot\-rnie.$cmName\.stk";
	    print "\tSTK output:       [$outFile]\n" if ($verbose);
	    my $modelString = '';
	    $modelString .= $modelsDir   if defined $modelsDir;
	    $modelString .= '/' . $model if defined $model;

	    system("cmalign -l -o $outFile $modelString " . $alnFileNames{$cmName} . '.fa') 
		and die "FATAL: failed to execute cmalign.\n[$!]";
	}
	
    }
    
}
else {
    die "FATAL: tabfile [$tabfile] is either missing or empty!\n";
}

exit(0);

######################################################################
#fileRoot: strip the suffix off file names 
sub fileRoot {
    my $file = shift;
    #print "file=[$file]\n" if (defined $verbose);
    #strip /'s from the filename:
    my @file = split(/\//, $file);
    $file = pop(@file) if (@file>1);
    #strip file suffix off:
    @file = split(/\./, $file);
    my $suff = pop(@file) if (scalar(@file)>1);
    $file = join('.', @file) if @file;
    #print "suffix=[$suff]\n" if (defined $verbose && defined $suff);
    #print "fileRoot=[$file]\n" if (defined $verbose);
    return $file;
}

######################################################################
#checkStrand: takes a star tand stop coord. Returns 1 if forward strand, -1 if reverse.
sub checkStrand{
    my ($a,$b)=@_;
    my $strand = 1;
    $strand = -1 if $b<$a;
    return $strand;
}

######################################################################
#overlapCheck: return true if the N/S-E overlaps with something in the %store hash:
sub overlapCheck{
    my ($name,$start,$end,$strand,$score,$evalue,$cmName,$store)=@_;
    return (0) if not defined $store->{$name};
    ($start,$end) = reorder($start,$end);
    my ($overlaps)=(0);
    for (my $i=0; $i<scalar(@{$store->{$name}}); $i++){
	my $a = $store->{$name}[$i]{'start'};
	my $b = $store->{$name}[$i]{'end'};
	my $s = checkStrand($a,$b);
	($a,$b) = reorder($a,$b);
	if (overlap($a,$b,$start,$end) && $strand == $s ){
	    my $overlapExtent = overlapExtent($a,$b,$start,$end);
	    if ($overlapExtent > 0.05){
		$overlaps=1;
		$store->{$name}[$i]{'cmNames'} .= ",$cmName=$score";
		#Replace entry if new score is higher:
		#Use e-value instead?
		if($store->{$name}[$i]{'bits'} < $score && $store->{$name}[$i]{'evalue'} > $evalue){
		    
		    $store->{$name}[$i]{'start'}    = $start;
		    $store->{$name}[$i]{'end'}      = $end;
		    $store->{$name}[$i]{'bits'}     = $score;
		    $store->{$name}[$i]{'evalue'}   = $evalue;
		    $store->{$name}[$i]{'sumBits'} += $score;
		    $store->{$name}[$i]{'cmName'}   = $cmName;
		}
		
	    }
	    
	}
    }
    return ($overlaps);
}

######################################################################
#Returns true if the coordinates for two regions ($x1, $y1) and ($x2, $y2) overlap:
# - assumes that $x1 < $y1 and $x2 < $y2.
sub overlap {
    my($x1, $y1, $x2, $y2) = @_;
    
    if ( ($x1<=$x2 && $x2<=$y1) || ($x1<=$y2 && $y2<=$y1) || ($x2<=$x1 && $x1<=$y2) || ($x2<=$y1 && $y1<=$y2)  ){
        return 1;
    }
    else {
        return 0;
    }
}

######################################################################
#Returns the extent of overlap between two regions A=($x1, $y1) and B=($x2, $y2):
#
# D = 2*|A n B|/(|A|+|B|)
#
sub overlapExtent {
    my($x1, $y1, $x2, $y2) = @_;
    
    ($x1, $y1) = reorder($x1, $y1);
    ($x2, $y2) = reorder($x2, $y2);
    # 1.
    # x1                   y1
    # |<---------A--------->|
    #    |<------B------>|
    #    x2             y2
    #    XXXXXXXXXXXXXXXXX
    # 2.  x1                     y1
    #     |<---------A----------->|
    # |<-------------B------>|
    # x2                    y2
    #     XXXXXXXXXXXXXXXXXXXX
    # 3. x1             y1
    #    |<------A------>|
    # |<---------B--------->|
    # x2                   y2
    #    XXXXXXXXXXXXXXXXX
    # 4. x1                    y1
    #    |<-------------A------>|
    #        |<---------B----------->|
    #        x2                     y2
    #        XXXXXXXXXXXXXXXXXXXX
    my $D=0;
    my $int=0;
    my $L1=$y1-$x1+1;
    my $L2=$y2-$x2+1;
    my $minL = min($L1,$L2);
    if ( ($x1<=$x2 && $x2<=$y1) && ($x1<=$y2 && $y2<=$y1) ){    #1.
	$D = $L2;
    }
    elsif ( ($x2<=$x1) && ($x1<=$y2 && $y2<=$y1) ){              #2.
	$D = $y2-$x1+1;
    }
    elsif ( ($x2<=$x1 && $x1<=$y2) && ($x2<=$y1 && $y1<=$y2) ){ #3.
	$D = $L1;
    }
    elsif ( ($x1<=$x2 && $x2<=$y1) && ($y1<=$y2) ){              #4.
	$D = $y1-$x2+1;
    }
    return $D/$minL;
}

######################################################################
#reorder: given 2 integers, return the smallest first & the largest last:
sub reorder {
    my ($x,$y)=@_;
    
    if ($y<$x){
	my $tmp = $x;
	$x = $y;
	$y = $tmp;
    }
    return ($x,$y);
}

######################################################################
#Max and Min
#max
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

#min
sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}

sub help {
    print STDERR <<EOF;

rnie.pl: 

Usage:   rnie.pl <options> -f genome.fasta
         rnie.pl <options> -t genome.tabfile

Options:       -h|--help                     Show this help.
               -v|--verbose                  Verbose mode. Prints lots of details to STDOUT.
	       
	       INPUT OPTIONS:
               -f|--fastafile <file>         Annotate Rho independent terminators on sequences in <file>.
               -t|--tabfile   <file>         Using a pre-computed tabfile <file> from Infernal. Resolves overlaps etc. 

	       RUNTIME OPTIONS:
	       -md|--modeldir <dir>          Directory containing the CM models [not required if the RNIE environment 
					     variable is set].
	       -m|--model     <cmfile>       Use cmfile rather than the default CM <\$RNIE/$modelGenome>
               -th|--thresh   <num>          Change the default bit-score threshold [Default=$cmsearchThreshDefaultGenome]
	       -e|--evalue    <num>          Use an E-value threshold rather than bit-scores

	       -o|--infOption <str>          Use additional infernal options. 
	       
	       RUN MODES:
	       --specific|--genome           Run in genome annotation mode. Uses a parameter set tuned for high 
	                                     specificity and speed [Default].
	       --sensitive|--gene            Run in targeted-region mode. Uses a parameter set tuned for high sensitivity. 
	                                     Primarily for answering the question: "Could there be a Rho independent 
                                             terminator following my gene of interest?"
					     
	       OUTPUT OPTIONS:
	       -p|--prefix <str>             Use 'str' as a prefix for the outputs rather than the default behaviour of 
	                                     stripping the suffix from fastafile/tabfile.

               -e|--embl                     Produce output in EMBL format
               -g|--gff                      Produce output in GFF format
	       -a|--aln                      Produce output in Stockholm (STK) format
	       
Dependencies:

Infernal version 1.0 or greater. The binaries for cmsearch, esl-sfetch
and cmalign must be in your path.

EXAMPLES:
To run the TRIT model:
rnie.pl -md <full path to model dir> -m trit.cm -f MTB_genome.fa 

Or if you have the \$RNIE environment variable set:
rnie.pl -m trit.cm -f MTB_genome.fa 

TODO: 

Add a Rho-dependent terminator model?:
1. find putative rut A and rut B sites using skew metrics (\#C-\#G)
2. Score spacing between cytosine residues -- find the 11\\pm1 pattern
3. Score the hairpin/boxB sequence
Can we use a profile HMM?
Is this pattern conserved across bacteria?

EOF
}




