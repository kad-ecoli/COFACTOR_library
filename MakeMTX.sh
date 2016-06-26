#!/bin/bash
DOCSTRING="blastpgp.sh protein.fasta [protein.blast [protein.mtx]]
   peform legacy psiblast on single FASTA sequence 'protein.fasta'
   output blast output to 'protein.blast' 
   output log-odd matrix to 'protein.mtx'"
if [ -z "$1" ];then
    echo "$DOCSTRING"
    exit
fi

#### parse input ####
blastdir="/home/zcx/Program/blast/2.2.26" # Directory for legacy BLAST executables
nrdb="/srv/zcx/nr/nr" # Non-redundant sequence database

fasta=`readlink -m "$1"`
prefix=`echo "$fasta"|sed "s/\.\w*$//g"`
if [ -z "$2" ];then
    blast_out="$prefix".blast
else
    blast_out=`readlink -m "$2"`
fi
if [ -z "$3" ];then
    mtx_out="$prefix".mtx
else
    mtx_out=`readlink -m "$3"`
fi

#### make working directory ####
tmpdir="/tmp/$USER/MakeMTX/$(basename $1)$(date '+%s')"
if [ ! -d $tmpdir ];then
    mkdir -p $tmpdir
else
    rm $tmpdir/*
fi
cd $tmpdir

#### perform psiblast ####
cp $fasta mtx.seq
$blastdir/bin/blastpgp -b 1000 -j 3 -h 0.001 -d $nrdb -i mtx.seq -C mtx.chk > blast.out
echo mtx.chk > mtx.pn
echo mtx.seq > mtx.sn
$blastdir/bin/makemat -P mtx  # make mtx.mtx

#### clean-up ####
cp blast.out "$blast_out"
cp mtx.mtx "$mtx_out"
rm -rf $tmpdir
