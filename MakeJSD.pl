#!/usr/bin/perl
use strict;
use File::Basename;
use Cwd 'abs_path';

my $docstring=qq{2011 by Ambrish Roy, 2016 by Chengxin Zhang
MakeJSD.pl protein.fasta blast.out JSD.dat SOP.dat
    protein.fasta - input FASTA sequence
    blast.out     - input legacy psi-blast output
    JSD.dat       - output Jensen–Shannon divergence
    SOP.dat       - output Sum of Pairs
};

if (! $ARGV[1] ){
    print $docstring;
    exit();
}

my $fasta    =abs_path($ARGV[0]);
my $blast_out=abs_path($ARGV[1]);
my $JSD_out  =abs_path($ARGV[2]);
my $SOP_out  =abs_path($ARGV[3]);
######################################################

## Scores each residues as putative binding site residue based on JSD and 
## PsiBlast output

my $s         =basename($fasta);
#my $datadir   ="!DATADIR!";  #for seq.txt and init.dat
my $user      =$ENV{LOGNAME} || $ENV{USER};
######################################################
my $tag          ="profile_$s"; 
my $workdir      ="/tmp/$user/$tag";
my $bindir  =dirname(abs_path(__FILE__));

`/bin/mkdir -p $workdir`;
chdir "$workdir";
`/bin/rm -fr $workdir/*`;

# CUSTOMIZE:

my $blastdir="/nfs/amino-library/blast"; # Directory containing  BLAST executables
my $nrdb    ="/nfs/amino-library/nr/nr";     # Non-redundant sequence database

#####################################################
# SETTINGS (Do not change below this line)
my @AA=
    (
     "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET","PHE",
     "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASX", "GLX", "UNK");

my %AA2index=
    (
     'A'=>'1', 'R'=>'2', 'N'=>'3', 'D'=>'4', 'C'=>'5', 'Q'=>'6', 'E'=>'7', 'G'=>'8', 'H'=>'9', 'I'=>'10', 'L'=>'11', 'K'=>'12', 
     'M'=>'13', 'F'=>'14', 'P'=>'15', 'S'=>'16', 'T'=>'17', 'W'=>'18', 'Y'=>'19', 'V'=>'20', 'B'=>'21', 'Z'=>'22', '-'=>'23');	
my %index2AA=
    (
     '1'=>'A','2'=>'R','3'=>'N','4'=>'D','5'=>'C','6'=>'Q','7'=>'E','8'=>'G','9'=>'H','10'=>'I','11'=>'L','12'=>'K', 
     '13'=>'M','14'=>'F','15'=>'P','16' =>'S','17'=>'T','18'=>'W','19'=>'Y','20'=>'V','21'=>'B','22'=>'Z','23'=>'-');	
	
	
     # A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
my @BLOSUM62 = 
    (
     [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0],
     [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3, 0, 1,-1],
     [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
     [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 1, 0,-1],
     [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2],
     [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
     [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 0, 2,-1],
     [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3, 0,-2,-1],
     [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 1, 0,-1],
     [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1],
     [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-3,-2,-1],
     [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1],
     [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-2, 0,-1],
     [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1],
     [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2],
     [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 1, 0, 0],
     [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0, 0,-1, 0],
     [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2],
     [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-2,-1,-1],
     [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1],
     [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 6, 0,-1],
     [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 5,-1],
     [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1]);

# This is the BLOSUM62 background distribution.
my %blsm_bck=
    (
     'A'=>0.078, 'R'=>0.051, 'N'=>0.041, 'D'=>0.052, 'C'=>0.024, 'Q'=>0.034, 'E'=>0.059, 'G'=>0.083, 
     'H'=>0.025, 'I'=>0.062, 'L'=>0.092, 'K'=>0.056, 'M'=>0.024, 'F'=>0.044, 'P'=>0.043, 'S'=>0.059, 
     'T'=>0.055, 'W'=>0.014, 'Y'=>0.034, 'V'=>0.072, 'B'=>0.041, 'Z'=>0.034, 'X'=>0.083,);
	 
my @BLOSM_FREQ = 
    ( [0.0215,0.0023,0.0019,0.0022,0.0016,0.0019,0.0030,0.0058,0.0011,0.0032,0.0044,0.0033,0.0013,0.0016,0.0022,0.0063,0.0037,0.0004,0.0013,0.0051,0.0019,0.0019,0.0058],
      [0.0023,0.0178,0.0020,0.0016,0.0004,0.0025,0.0027,0.0017,0.0012,0.0012,0.0024,0.0062,0.0008,0.0009,0.0010,0.0023,0.0018,0.0003,0.0009,0.0016,0.0020,0.0025,0.0017],
      [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
      [0.0022,0.0016,0.0037,0.0213,0.0004,0.0016,0.0049,0.0025,0.0010,0.0012,0.0015,0.0024,0.0005,0.0008,0.0012,0.0028,0.0019,0.0002,0.0006,0.0013,0.0037,0.0016,0.0025],
      [0.0016,0.0004,0.0004,0.0004,0.0119,0.0003,0.0004,0.0008,0.0002,0.0011,0.0016,0.0005,0.0004,0.0005,0.0004,0.0010,0.0009,0.0001,0.0003,0.0014,0.0004,0.0003,0.0008],
      [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
      [0.0030,0.0027,0.0022,0.0049,0.0004,0.0035,0.0161,0.0019,0.0014,0.0012,0.0020,0.0041,0.0007,0.0009,0.0014,0.0030,0.0020,0.0003,0.0009,0.0017,0.0022,0.0035,0.0019],
      [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378],
      [0.0011,0.0012,0.0014,0.0010,0.0002,0.0010,0.0014,0.0010,0.0093,0.0006,0.0010,0.0012,0.0004,0.0008,0.0005,0.0011,0.0007,0.0002,0.0015,0.0006,0.0014,0.0010,0.0010],
      [0.0032,0.0012,0.0010,0.0012,0.0011,0.0009,0.0012,0.0014,0.0006,0.0184,0.0114,0.0016,0.0025,0.0030,0.0010,0.0017,0.0027,0.0004,0.0014,0.0120,0.0010,0.0009,0.0014],
      [0.0044,0.0024,0.0014,0.0015,0.0016,0.0016,0.0020,0.0021,0.0010,0.0114,0.0371,0.0025,0.0049,0.0054,0.0014,0.0024,0.0033,0.0007,0.0022,0.0095,0.0014,0.0016,0.0021],
      [0.0033,0.0062,0.0024,0.0024,0.0005,0.0031,0.0041,0.0025,0.0012,0.0016,0.0025,0.0161,0.0009,0.0009,0.0016,0.0031,0.0023,0.0003,0.0010,0.0019,0.0024,0.0031,0.0025],
      [0.0013,0.0008,0.0005,0.0005,0.0004,0.0007,0.0007,0.0007,0.0004,0.0025,0.0049,0.0009,0.0040,0.0012,0.0004,0.0009,0.0010,0.0002,0.0006,0.0023,0.0005,0.0007,0.0007],
      [0.0016,0.0009,0.0008,0.0008,0.0005,0.0005,0.0009,0.0012,0.0008,0.0030,0.0054,0.0009,0.0012,0.0183,0.0005,0.0012,0.0012,0.0008,0.0042,0.0026,0.0008,0.0005,0.0012],
      [0.0022,0.0010,0.0009,0.0012,0.0004,0.0008,0.0014,0.0014,0.0005,0.0010,0.0014,0.0016,0.0004,0.0005,0.0191,0.0017,0.0014,0.0001,0.0005,0.0012,0.0009,0.0008,0.0014],
      [0.0063,0.0023,0.0031,0.0028,0.0010,0.0019,0.0030,0.0038,0.0011,0.0017,0.0024,0.0031,0.0009,0.0012,0.0017,0.0126,0.0047,0.0003,0.0010,0.0024,0.0031,0.0019,0.0038],
      [0.0037,0.0018,0.0022,0.0019,0.0009,0.0014,0.0020,0.0022,0.0007,0.0027,0.0033,0.0023,0.0010,0.0012,0.0014,0.0047,0.0125,0.0003,0.0009,0.0036,0.0022,0.0014,0.0022],
      [0.0004,0.0003,0.0002,0.0002,0.0001,0.0002,0.0003,0.0004,0.0002,0.0004,0.0007,0.0003,0.0002,0.0008,0.0001,0.0003,0.0003,0.0065,0.0009,0.0004,0.0002,0.0002,0.0004],
      [0.0013,0.0009,0.0007,0.0006,0.0003,0.0007,0.0009,0.0008,0.0015,0.0014,0.0022,0.0010,0.0006,0.0042,0.0005,0.0010,0.0009,0.0009,0.0102,0.0015,0.0007,0.0007,0.0008],
      [0.0051,0.0016,0.0012,0.0013,0.0014,0.0012,0.0017,0.0018,0.0006,0.0120,0.0095,0.0019,0.0023,0.0026,0.0012,0.0024,0.0036,0.0004,0.0015,0.0196,0.0012,0.0012,0.0018],
      [0.0019,0.0020,0.0141,0.0037,0.0004,0.0015,0.0022,0.0029,0.0014,0.0010,0.0014,0.0024,0.0005,0.0008,0.0009,0.0031,0.0022,0.0002,0.0007,0.0012,0.0141,0.0015,0.0029],
      [0.0019,0.0025,0.0015,0.0016,0.0003,0.0073,0.0035,0.0014,0.0010,0.0009,0.0016,0.0031,0.0007,0.0005,0.0008,0.0019,0.0014,0.0002,0.0007,0.0012,0.0015,0.0073,0.0014],
      [0.0058,0.0017,0.0029,0.0025,0.0008,0.0014,0.0019,0.0378,0.0010,0.0014,0.0021,0.0025,0.0007,0.0012,0.0014,0.0038,0.0022,0.0004,0.0008,0.0018,0.0029,0.0014,0.0378]);



#####################################################
my %AA2sim=();
for(my $i=0;$i<23;$i++)
{
    for(my $j=0;$j<23;$j++)
    {
	my $x= $i+1;my $y= $j+1;
	$AA2sim{$index2AA{$x},$index2AA{$y}}= $BLOSUM62[$i][$j];
    }
}
my $PSEUDOCOUNT= 0.0000001;
my $lambda     = 0.5;
my $window     = 3;

###################################################################################################
#`/bin/cp $bindir/alignhits.pl   .`;

my @seqtxts=`cat $fasta`;
my $sequence="";    
my %seqQ=();
foreach my $seqtxt(@seqtxts)
{
	goto pos1 if($seqtxt=~/\>/);
	$seqtxt=~s/\s//mg;
	$seqtxt=~s/\n//mg;
	$sequence=$sequence.$seqtxt;
  pos1:;
}
my $Lch=length $sequence;  

open(FASTA,">protein.fasta");
printf FASTA ">protein\n";
for(my $i=1;$i<=$Lch;$i++)
{
	my $seq1=substr($sequence,$i-1,1);
	$seqQ{$i}=$seq1;
	printf FASTA "$seq1";
	if(int($i/60)*60==$i)
	{
		printf FASTA "\n";
	}
}
close(FASTA);

if(!-s "$blast_out")
{
	system ("$blastdir/bin/blastpgp  -b 1000 -j 3 -h 0.001 -d $nrdb -i protein.fasta > blast.out");
	`/bin/cp  blast.out $blast_out`;	
}
else
{
	`/bin/cp  $blast_out  blast.out `;
}

`$bindir/alignhits.pl -fas blast.out -q protein.fasta blast.fasta`;

my %ALN=();my $Pcount=0;
open(ALN,"<$workdir/blast.fasta") || die "Cant open blast.fasta";
while(my $line=<ALN>)
{
    chomp($line);
    if($line=~/^>(\S+)/)
    {
	my $Pname=$1;
	$Pcount++;
	my $Evalue= $1 if($line=~/E=(\S+)/);	
	$ALN{$Pcount,0}=$Pname;
	$ALN{$Pcount,1}=$Evalue;
    }
    else
    {
	$line=~s/X/-/g;
	$ALN{$Pcount,2}=$line;	
    }    
}
close(ALN);

if($Pcount > 10)
{
   $Pcount=1000 if ($Pcount > 1000);
   my %seq_weights = ();
   &load_sequence_weights(\%ALN,\$Pcount,\%AA2index,\%seq_weights);

   my %jsd_scores = ();my %zscore1=(); my %wscore1=();
   &js_divergence_score(\%ALN,\$Pcount,\%blsm_bck,\%AA2index,\%seq_weights,$PSEUDOCOUNT,$lambda,\%jsd_scores);
   &z_scores(\%jsd_scores,\$Lch,\%zscore1);

   #### Window removed because of low performance
   #&w_scores(\%jsd_scores,\$Lch,\$window,\$lambda,\%wscore1);
   #&z_scores(\%wscore1,\$Lch,\%zscore1);

   my %SOP=();my %zscore2=();my %wscore2=();
   &sum_of_pairs(\%ALN,\$Pcount,\%AA2sim,\%seq_weights,\%SOP);
   &z_scores(\%SOP,\$Lch,\%zscore2);

   #### Window removed because of low performance
   #&w_scores(\%SOP,\$Lch,\$window,\$lambda,\%wscore2);
   #&z_scores(\%wscore2,\$Lch,\%zscore2);

   open(OUT1,">$workdir/JSD.dat");
   open(OUT2,">$workdir/SOPM.dat");
   for(my $j=1;$j<=$Lch;$j++)
   {
     printf OUT1 "%-4d  %5s\n",$j,$zscore1{$j};
     printf OUT2 "%-4d  %5s\n",$j,$zscore2{$j};    
   }
   close(OUT1);
   close(OUT2);
}
else
{
   open(OUT1,">$workdir/JSD.dat");
   open(OUT2,">$workdir/SOPM.dat");
   for(my $j=1;$j<=$Lch;$j++)
   {
     printf OUT1 "%-4d  %5s\n",$j,'1.00';
     printf OUT2 "%-4d  %5s\n",$j,'1.00';
   }
   close(OUT1);
   close(OUT2);
}
`/bin/cp $workdir/JSD.dat  $JSD_out`;
`/bin/cp $workdir/SOPM.dat $SOP_out`;

`sync`;
my $time=`date`;
printf "Ended: $time";
sleep(1);
`rm -fr $workdir`;


sub sum_of_pairs
{
    my ($ALN_ref,$Nseq_ref,$blossum_ref,$SW_ref,$SOP)=@_;
    my %align       =%$ALN_ref;
    my $Nseq        =$$Nseq_ref;
    my %sim_matrix  =%$blossum_ref;
    my %weights     =%$SW_ref;  
   
    my @Qres        =split(//,$align{1,2});
    my $Ncol        =$#Qres;    
    my $Qresno=0;my %Qmapping=();
    for(my $j=0;$j<=$#Qres;$j++)
    {
	if($Qres[$j] ne '-'){$Qresno++;$Qmapping{$Qresno}=$j;} 	
    }
    my @ARR=();my $sum_seq_weights = 0;
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	for(my $j=0;$j<=$#res;$j++)
	{
	    $ARR[$i][$j]=$res[$j];
	}
	$sum_seq_weights += $weights{$i};
    }
    my %SOP_scores=();
    for(my $C=0;$C<=$Ncol;$C++)
    {
	my $sum   = 0;
	my $maxsum= 0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	    my $aa1=$ARR[$i][$C];
	    for(my $j=$i+1;$j<=$Nseq;$j++)
	    {
		my $aa2=$ARR[$j][$C];
		if(($ARR[$i][$C] eq '-')||($ARR[$j][$C] eq '-')){next;}
		$maxsum += $weights{$i} * $weights{$j};
		$sum    += $weights{$i} * $weights{$j} * $sim_matrix{$aa1,$aa2};	
	    }	    
	}
	if ($maxsum !=0)
	{
	    $sum   /=$maxsum;
	}
	else
	{
	    $sum=0;
	}
	my $gap_sum=0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	  if($ARR[$i][$C] eq '-')
	  {
	      $gap_sum +=   $weights{$i}; 
	  }
	}
	my $weighted_gap_penalty = (1-  ($gap_sum/$sum_seq_weights));
	$SOP_scores{$C}=($sum*$weighted_gap_penalty);	
    }
    for(my $j=1;$j<=$Qresno;$j++)
    {
	$$SOP{$j}= sprintf("%5.4f",$SOP_scores{$Qmapping{$j}});
    }
    return;
    
}
sub js_divergence_score
{
    my ($ALN_ref,$Nseq_ref,$dist_ref,$AA_ref,$SW_ref,$pseudocount,$lambda,$JSD_ref)=@_;
    my %align  =%$ALN_ref;
    my $Nseq   =$$Nseq_ref;
    my %distr  =%$dist_ref;
    my %weights=%$SW_ref;  
    my %AA2in  =%$AA_ref;    
    my @Qres   =split(//,$align{1,2});
    my $Ncol   =$#Qres;
    
    my $Qresno=0;my %Qmapping=();
    for(my $j=0;$j<=$#Qres;$j++)
    {
	if($Qres[$j] ne '-'){$Qresno++;$Qmapping{$Qresno}=$j;} 	
    }  
    my @ARR=();
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	for(my $j=0;$j<=$#res;$j++)
	{
	    $ARR[$i][$j]=$res[$j];
	}
    }
    my $AAcount = keys %AA2in;
    my %AA_freq=();my %sum_seq_weights=();my %R=();my %JSD_scores=();
    for(my $j=0;$j<=$Ncol;$j++)
    {
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{
	    $AA_freq{$j,$key}=0;
	    $R{$j,$key}      =0;
	}
	$sum_seq_weights{$j} =0;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	    $AA_freq{$j,$ARR[$i][$j]} += (1.0 * $weights{$i});
	    $sum_seq_weights{$j}      += $weights{$i}; 
	}
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{	   
	    $AA_freq{$j,$key}=$AA_freq{$j,$key}/($sum_seq_weights{$j} + $AAcount * $pseudocount);
	    if($key ne '-')
	    {
		$R{$j,$key}=($lambda*$AA_freq{$j,$key}) + ((1-$lambda)*$distr{$key});	   
	    }
	}
    }
    for(my $j=0;$j<=$Ncol;$j++)
    {
	my $JSD_score=0;my $weighted_gap_penalty=0;my $gap_sum=0;
	foreach my $key (sort {$AA2in{$a} <=> $AA2in{$b}} keys %AA2in)
	{
	    if($key eq '-'){next;}
	    if($R{$j,$key} != 0.0)
	    {
		if($AA_freq{$j,$key} == 0.0)
		{
		    $JSD_score += ($distr{$key} * log2($distr{$key}/$R{$j,$key}));
		}
		elsif($distr{$key} == 0.0)
		{
		    $JSD_score += ($AA_freq{$key} * log2($AA_freq{$key}/$R{$j,$key}));
		}
		else
		{
		    $JSD_score += ($AA_freq{$j,$key} * log2($AA_freq{$j,$key}/$R{$j,$key})) + ($distr{$key} * log2($distr{$key}/$R{$j,$key}));
		}
	    }
	}	
	$JSD_score = $JSD_score/2;
	for(my $i=1;$i<=$Nseq;$i++)
	{
	  if($ARR[$i][$j] eq '-')
	  {
	      $gap_sum +=   $weights{$i}; 
	  }
	}
	$weighted_gap_penalty = (1-  ($gap_sum/$sum_seq_weights{$j}));
	$JSD_scores{$j}=($JSD_score*$weighted_gap_penalty);	
    }
    for(my $j=1;$j<=$Qresno;$j++)
    {
	$$JSD_ref{$j}= sprintf("%5.4f",$JSD_scores{$Qmapping{$j}});
    }
    return;
}

# Get Seqeuence Weights based on Henikoff-Henikoff schema
sub load_sequence_weights
{
    my ($ALN_ref,$Nseq_ref,$AA_ref,$seq_weights_ref)=@_;

    my %align  =%$ALN_ref;
    my $Nseq   =$$Nseq_ref;
    my %AA2in  =%$AA_ref;

    my %NRC=();my %RC=();my %seen=();
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});

	for(my $j=0;$j<=$#res;$j++)
	{
	    my $AAN=$AA2in{$res[$j]};
	    $RC{$j,$AAN}++;	  
	    if(exists $seen{$j,$AAN}){next;}
	    else{$seen{$j,$AAN}=1;$NRC{$j}++;}	    
	}	
    }   
    for(my $i=1;$i<=$Nseq;$i++)
    {
	my @res=split(//,$align{$i,2});
	
	for(my $j=0;$j<=$#res;$j++)
	{
	    my $AAN=$AA2in{$res[$j]};
	    $$seq_weights_ref{$i} += 1/($NRC{$j} * $RC{$j,$AAN});    	   
	}
	$$seq_weights_ref{$i} = sprintf("%6.4f",$$seq_weights_ref{$i}/$#res);
    }
    return;  
}
	
sub log2 
{
    my $n = shift;
    return log($n)/log(2);
}

sub w_scores
{# This subroutine takes a list of scores and a length and transforms them so that each position is a weighted average of the surrounding positions. 
     my($score_ref,$L,$W,$lambda,$w_score_ref)=@_;
     my %r_score  =%$score_ref;
     my $Lch   =$$L;
     my $window   =$$W;
     for(my $i=1;$i<=$Lch;$i++)
     {
	 my $sum     = 0;
	 my $n_terms = 0;
	 for(my $j=$i-$window;$j<=$i+$window;$j++)
	 {
	     if(($i != $j) && (exists $r_score{$j}))
	     {
		 $sum +=  $r_score{$j};
		 $n_terms++;
	     }
	 }
	 $$w_score_ref{$i}= sprintf("%5.4f",$sum/$n_terms);
     }
}
    
sub z_scores
{# This subroutine calculates Z-scores for a set of scores (NOTE: hash index starts from 1)   
     my($score_ref,$L,$z_score_ref)=@_;
     my %r_score  =%$score_ref;
     my $Lch   =$$L;
     my $z1_a  = 0; my $z1_a2=0; my %zscore=();
     for(my $i=1;$i<=$Lch;$i++)
     {	 
	 $z1_a  += $r_score{$i}; 
	 $z1_a2 += $r_score{$i}**2;
     }
     $z1_a/=$Lch;
     $z1_a2/=$Lch;
     my $z1_sd=sqrt($z1_a2-$z1_a**2);
     for(my $i=1;$i<=$Lch;$i++)
     {
         $$z_score_ref{$i}=($r_score{$i}-$z1_a)/$z1_sd;
     }
}

