#!/usr/bin/perl
use strict;
my $docstring=qq{MakePRF.pl protein.fasta blast.out > protein.prf
    calculate Henikoff weight "protein.prf" according to FASTA input 
    sequence "protein.fasta" & legacy psi-blast output "blast.out"
};

if ( ! $ARGV[1] ){
    print $docstring;
    exit();
}

my $fasta=$ARGV[0];
my $blast_out=$ARGV[1];
my $prf=$ARGV[2];

my @AA=('-','A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y');

my @seqtxts=`cat $fasta`;
my $sequence="";
my %seqQ=();
foreach my $seqtxt(@seqtxts){
    if($seqtxt!~/\>/){
        $seqtxt=~s/\s+//mg;
        $seqtxt=~s/\n//mg;
        $sequence=$sequence.$seqtxt;
    }
}
my $Lch=length $sequence;

for(my $i=1;$i<=$Lch;$i++){
    my $a=substr($sequence,$i-1,1);
    $seqQ{$i}=$a;   #only for check
}

my $ROUND=1;my $it=0;my %am=();my %nA=();my %w=();my $w_all=0; my %ev_dim=();
open(BL,"$blast_out");
while(my $line=<BL>)
{
    if($line=~/Results from round\s+(\d+)/){
	$ROUND=$1;
    }
}
seek(BL,0,0);
while(my $line=<BL>)
{
    if($line=~/round\s+$ROUND/){
	while(my $line2=<BL>){
	    if($line2=~/Expect =\s*(\S+)/){
		my $ev=$1;
		<BL>=~/Identities =\s*\S+\s+\((\S+)\%/;
		<BL>;
		my $id=$1;
		if($ev<0.001){
		    $it++; #number of aligned sequences
		    $ev_dim{$it}=$ev;
		    while(my $line3=<BL>){
			if($line3=~/Query\:\s*(\d+)\s+(\S+)\s+(\S+)/){
			    my $i1=$1;
			    my $seq1=$2;
			    <BL>;
			    <BL>=~/Sbjct\:\s*(\S+)\s+(\S+)\s+(\S+)/;
			    my $seq2=$2;
			    <BL>;
			    ###
			    my $L=length $seq1;
			    my $m1=$i1-1;
			    for(my $i=1;$i<=$L;$i++){
				my $q1=substr($seq1,$i-1,1);
				my $q2=substr($seq2,$i-1,1);
				$m1++ if($q1 ne '-');
				$am{$it,$m1}=$q2 if($q1 ne '-' && $q2 ne '-');
			    }
			}
			else{
			    last;#goto pos2;   
			}
		    }
		}
	      pos2:;
	    }
	}
    }
}		
close(BL);
$it++;
for(my $i=1;$i<=$Lch;$i++){
    $am{$it,$i}=$seqQ{$i};
}
####### Henikoff-Henikoff weight $w{i_seq} ----------->
##### nA{A,i_pos}: number of times A appear at the position:

for(my $i=1;$i<=$Lch;$i++){
    for(my $j=1;$j<=$it;$j++){
	$nA{$am{$j,$i},$i}++;
    }
}
##### henikoff weight w(i)=sum of 1/rs:
for(my $i=1;$i<=$it;$i++)
{
    for(my $j=1;$j<=$Lch;$j++){
	####### r: number of different residues at j'th position:
	my $r=0;
	foreach my $A(@AA){
	    $r++ if($nA{$A,$j}>0);
	}
	my $A=$am{$i,$j};
	my $s2=$nA{$A,$j};
	my $D=($r*$s2);
	if($D > 0)
	{
		my $w1=1.0/($r*$s2);
		$w{$i}+=$w1;
	}
    }
    $w_all+=$w{$i};
}
#### normalization of w(i):
for(my $i=1;$i<=$it;$i++){
    $w{$i}/=$w_all;
}
#^^^^^ Henikoff weight finished ^^^^^^^^^^^^^^^^^
###########Henikoff  weighted frequency #################
my %log=();;
for(my $i=1;$i<=$it;$i++)
{
    for(my $j=1;$j<=$Lch;$j++)
    {
	my $A=$am{$i,$j};
	$log{$j,$A}+=$w{$i};
    }
}

my $txt="$Lch\n";
for(my $j=1;$j<=$Lch;$j++){
    $txt.=sprintf("%3d $seqQ{$j} %3d",$j,$j);
    my $norm=0;
    for(my $i=1;$i<=22;$i++){
	$norm +=$log{$j,$AA[$i]};
    }
    for(my $i=1;$i<=22;$i++){
	$txt.=sprintf("%10.3f",$log{$j,$AA[$i]}/$norm);
    }
    $txt.="\n";
}

if ($prf){
    open(FREQ,">$prf");
    print FREQ "$txt";
    close(FREQ);
}else{
    print "$txt";
}

exit();
