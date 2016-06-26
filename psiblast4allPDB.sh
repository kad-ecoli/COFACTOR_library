#!/bin/bash
DOCSTRING="psiblast4allPDB.sh PDBDIR/ list
   Given 'PDBDIR', a folder of single chain PDB files, 
   and 'list', a list of PDB chains. Seqentially run
   fixMSE.py reindex_pdb.py getfasta.py MakeMTX.sh MakePRF.pl MakeJSD.pl
   Deposit *.fasta files at PDBDIR/FASTA/
   Deposit *.mtx files at PDBDIR/MTX/
   Deposit *.prf files at PDBDIR/PRF/
   Deposit *.blast, *sop, *jsd files at PDBDIR/JSD/

psiblast4allPDB.sh -runstyle=parallel PDBDIR/ list
   qsub jobs for each chain

psiblast4allPDB.sh -runstyle=gnuparallel PDBDIR/ list
   use gnu parallel to parallelize jobs for each chain

psiblast4allPDB.sh -outfmt=sequence PDBDIR/ list
   only generate the FASTA format sequence database
"
if [ -z "$2" ];then
    echo "$DOCSTRING"
    exit
fi

#### parse input ####
RUNSTYLE="serial"
OUTFMT="cofactor"
if [ "$1" == "-runstyle=parallel"  ];then
    RUNSTYLE="parallel"
    shift
elif [ "$1" == "-runstyle=gnuparallel"  ];then
    RUNSTYLE="gnuparallel"
    shift
fi
if [ "$1" == "-outfmt=sequence"  ];then
    OUTFMT="sequence"
    shift
fi
PDBDIR="$(readlink -m $1)"
LIST="$(readlink -m $2)"
BINDIR=`dirname "$(readlink -m $0)"`

if [ ! -d "$PDBDIR" ];then
    echo "ERROR! No such directory: $PDBDIR"
    exit
elif [ ! -f "$LIST" ];then
    echo "ERROR! No such file: $LIST"
    exit
elif [ "$OUTFMT" == "sequence" ];then
    FASTADIR="$(readlink -m $3)"
    for chain in `cat $LIST|sed 's/\.pdb//g'`;do
        if [ -s "$PDBDIR/$chain.pdb" ];then
	    if [ -z "$FASTADIR" ];then
                $BINDIR/pdb2fasta.py -outfmt=COFACTOR $PDBDIR/$chain.pdb
	    elif [ ! -s "$FASTADIR/$chain.fasta" ];then
                $BINDIR/pdb2fasta.py -outfmt=COFACTOR $PDBDIR/$chain.pdb > $FASTADIR/$chain.fasta
	    fi
        fi
    done
    exit
fi

#### parallel run by PBS ####
WALLTIME="mem=4000m,walltime=4:00:00"
Q="default" # queue destination
#Q="urgent" # queue destination
N_JOB_MAX="300" # maximum number of jobs submitted by user
N_JOB_TOTAL_MAX="1200" # maximum number of jobs submitted by everyone

if [ $RUNSTYLE == "parallel" ];then
    mkdir -p record/
    RECORDDIR="$(readlink -m record/)"
    for chain in `cat $LIST|sed 's/\.pdb//g'`;do
        printf $chain

	## check if job finished ##
        if [ -s "$PDBDIR/FASTA/$chain.fasta" ] && \
	   [ -s "$PDBDIR/MTX/$chain.mtx" ]     && \
	   [ -s "$PDBDIR/PRF/$chain.prf" ]     && \
	   [ -s "$PDBDIR/JSD/$chain.jsd" ]     && \
	   [ -s "$PDBDIR/JSD/$chain.sop" ];then
	    printf " completed\n"
	    continue
	fi
	printf "\n"

        ## write PBS script ##
	echo "$chain" > $RECORDDIR/$chain
	JOBNAME="CFl_$chain"
        cat > "$RECORDDIR/$JOBNAME.pbs" <<EOF
#!/bin/bash
#PBS -N $JOBNAME
#PBS -l $WALLTIME
#PBS -q $Q
#PBS -e $RECORDDIR/$JOBNAME.err
#PBS -o $RECORDDIR/$JOBNAME.out
#PBS -d $PDBDIR
$(readlink -m $0) $PDBDIR $RECORDDIR/$chain
EOF
        chmod +x "$RECORDDIR/$JOBNAME.pbs"

	## qsub PBS script ##
        QSUB=""
	until [ ! -z "$QSUB" ];do
	    N_JOB_RUN="$(qstat|grep $USER|wc -l)"
	    N_JOB_RUN_TOTOAL="$(qstat |wc -l)"
	    if [ ! -z "$(qstat|grep $JOBNAME)" ];then
	        QSUB="running"
	    elif [ $N_JOB_RUN -lt $N_JOB_MAX ] || [ $N_JOB_RUN_TOTOAL -lt $N_JOB_TOTAL_MAX ];then
	        QSUB=`qsub "$RECORDDIR/$JOBNAME.pbs"`
	    else
	        printf '.'
	        sleep 30
	    fi
	done
        printf " $QSUB\n"
    done
    exit
fi

#### parallel run by gnu parallel ####
CPU_NUM=20 # number of cores to use
if [ $RUNSTYLE == "gnuparallel" ];then
    mkdir -p record/
    RECORDDIR="$(readlink -m record/)"
    JOBNAME="gnuparallel.sh"
    cat > "$RECORDDIR/$JOBNAME" <<EOF
#!/bin/bash
read -r -d '' CMD_LST << EOF
EOF
    chmod +x "$RECORDDIR/$JOBNAME"
    for chain in `cat $LIST|sed 's/\.pdb//g'`;do
        printf $chain

	## check if job finished ##
        if [ -s "$PDBDIR/FASTA/$chain.fasta" ] && \
	   [ -s "$PDBDIR/MTX/$chain.mtx" ]     && \
	   [ -s "$PDBDIR/PRF/$chain.prf" ]     && \
	   [ -s "$PDBDIR/JSD/$chain.jsd" ]     && \
	   [ -s "$PDBDIR/JSD/$chain.sop" ];then
	    printf " completed\n"
	    continue
	fi
	printf "\n"

        ## write PBS script ##
	echo "$chain" > $RECORDDIR/$chain
        echo "$(readlink -m $0) $PDBDIR $RECORDDIR/$chain" >> "$RECORDDIR/$JOBNAME"
    done
    echo  "EOF" >> "$RECORDDIR/$JOBNAME"
    echo  "echo \"\$CMD_LST\"|time parallel -k -j $CPU_NUM" >> "$RECORDDIR/$JOBNAME"
    `$RECORDDIR/$JOBNAME`
    exit
fi

#### make directory ####
mkdir -p "$PDBDIR/FASTA"
mkdir -p "$PDBDIR/MTX"
mkdir -p "$PDBDIR/PRF"
mkdir -p "$PDBDIR/JSD"
#mkdir -p "$PDBDIR/sequence"

TMPDIR="/tmp/$USER/psiblast4allPDB_$(basename $LIST)"
if [ -d "$TMPDIR" ];then
    rm -rf $TMPDIR/*
fi

mkdir -p "$TMPDIR/FASTA"
mkdir -p "$TMPDIR/MTX"
mkdir -p "$TMPDIR/PRF"
mkdir -p "$TMPDIR/JSD"
mkdir -p "$TMPDIR/bin"
#mkdir -p "$TMPDIR/sequence"

#### copy executables ####
fixMSE="$TMPDIR/bin/fixMSE.py"
reindex_pdb="$TMPDIR/bin/reindex_pdb.py"
pdb2fasta="$TMPDIR/bin/pdb2fasta.py"
alignhits="$TMPDIR/bin/alignhits.pl"
MakeMTX="$TMPDIR/bin/MakeMTX.sh"
MakePRF="$TMPDIR/bin/MakePRF.pl"
MakeJSD="$TMPDIR/bin/MakeJSD.pl"

cp "$BINDIR/fixMSE.py"      $fixMSE
cp "$BINDIR/reindex_pdb.py" $reindex_pdb
cp "$BINDIR/pdb2fasta.py"   $pdb2fasta
cp "$BINDIR/alignhits.pl"   $alignhits
cp "$BINDIR/MakeMTX.sh"     $MakeMTX
cp "$BINDIR/MakePRF.pl"     $MakePRF
cp "$BINDIR/MakeJSD.pl"     $MakeJSD

chmod +x $TMPDIR/bin/*

#### fixMSE pdb2fasta MakeMTX MakePRF MakeJSD ####
for chain in `cat $LIST|sed 's/\.pdb//g'`;do
    if [ -s "$PDBDIR/FASTA/$chain.fasta" ] && \
       [ -s "$PDBDIR/MTX/$chain.mtx" ]     && \
       [ -s "$PDBDIR/PRF/$chain.prf" ]     && \
       [ -s "$PDBDIR/JSD/$chain.jsd" ]     && \
       [ -s "$PDBDIR/JSD/$chain.sop" ];then
        continue # every data file is generated for this target
    fi

    cp $PDBDIR/$chain.pdb $TMPDIR/$chain.pdb
    printf $chain

    # fixMSE
    $fixMSE -clean=true $TMPDIR/$chain.pdb $TMPDIR/$chain.pdb
    cp $TMPDIR/$chain.pdb $PDBDIR/$chain.pdb
    printf " fixMSE"

    # reindex_pdb
    $reindex_pdb 1 $TMPDIR/$chain.pdb $TMPDIR/$chain.pdb
    cp $TMPDIR/$chain.pdb $PDBDIR/$chain.pdb
    printf " reindex_pdb"

    # pdb2fasta
    $pdb2fasta -outfmt=COFACTOR $TMPDIR/$chain.pdb > $TMPDIR/FASTA/$chain.fasta
    if [ ! -s "$PDBDIR/FASTA/$chain.fasta" ];then
        cp "$TMPDIR/FASTA/$chain.fasta" "$PDBDIR/FASTA/$chain.fasta"
        printf " pdb2fasta"
    fi

    # MakeMTX
    if [ ! -s "$PDBDIR/MTX/$chain.mtx" ]||[ ! -s "$PDBDIR/JSD/$chain.blast" ];then
	$MakeMTX $TMPDIR/FASTA/$chain.fasta $TMPDIR/JSD/$chain.blast $TMPDIR/MTX/$chain.mtx
	cp "$TMPDIR/MTX/$chain.mtx"   "$PDBDIR/MTX/$chain.mtx"
	#cp "$TMPDIR/JSD/$chain.blast" "$PDBDIR/JSD/$chain.blast" # do not copy blast file back to save space
	printf " MakeMTX"
    else
	cp "$PDBDIR/JSD/$chain.blast" "$TMPDIR/JSD/$chain.blast"
    fi

    # MakePRF
    if [ ! -s "$PDBDIR/PRF/$chain.prf" ];then
	$MakePRF $TMPDIR/FASTA/$chain.fasta $TMPDIR/JSD/$chain.blast $TMPDIR/PRF/$chain.prf
	cp "$TMPDIR/PRF/$chain.prf" "$PDBDIR/PRF/$chain.prf"
	printf " MakePRF"
    fi

    # MakeJSD
    if [ ! -s "$PDBDIR/JSD/$chain.jsd" ]||[ ! -s "$PDBDIR/JSD/$chain.sop" ];then
	$MakeJSD $TMPDIR/FASTA/$chain.fasta $TMPDIR/JSD/$chain.blast $TMPDIR/JSD/$chain.jsd $TMPDIR/JSD/$chain.sop
	cp "$TMPDIR/JSD/$chain.jsd" "$PDBDIR/JSD/$chain.jsd"
	cp "$TMPDIR/JSD/$chain.sop" "$PDBDIR/JSD/$chain.sop"
	printf " MakeJSD"
    fi

    printf "\n"
done

#### cat all fasta files ####
#cat $TMPDIR/FASTA/*.fasta > $TMPDIR/sequence/lib.fasta
#cp "$TMPDIR/sequence/lib.fasta" "$PDBDIR/sequence/lib.fasta"

#### clean up ####
rm -rf $TMPDIR
