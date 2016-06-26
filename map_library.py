#!/usr/bin/env python
docstring='''
map_library.py A.fasta B.fasta > mapAB
    blast using each entry in B.fasta as query and A.fasta as database.
    Output "mapAB", a mapping between entries in "A.fasta" and "B.fasta"
    Each line consists of three fields:
    "one_entry_in_A" "num_match_in_B" "list_of_matched_B_entry"

options:
-ncbidir=/nfs/amino-library/blast
    search formatdb and blastall under ncbidir/bin

-phmmer_path=phmmer
    use phmmer at phmmer_path instead of blastall for mapping

-mapfirst={true,false} how many blastp hits to map
    true  - (default) only map the first blastp hit
    false - map all hits
'''

import os,sys
import shutil
import subprocess
import re
import random # for creating random temporary file

hit_pattern=re.compile("\n(\w+)")
hit_txt="Sequences producing significant alignments"

def split_fasta(infile):
    '''Split multi-sequence FASTA file "infile" into a list for all single-
    sequence fasta-format text. Also return a list of all headers. 
    Each header only preserve the first word'''
    fp=open(infile,'rU')
    txt=fp.read()
    fp.close()
    fasta_list=['>'+'\n'.join([line.split()[0] for line in e.splitlines()]
        )+'\n' for e in txt.split('>') if e.strip()]
    header_list=[e.splitlines()[0][1:] for e in fasta_list]
    return fasta_list,header_list

def map_library(fastaA,fastaB,ncbidir='',mapfirst=True):
    '''create blast mapping between two fasta sequence database "fastaA" 
    and "fastaB". Return mapping as text.
    "ncbidir" is the path to formatdb and to blastall
    "mapfirst" whether only take the first blastp hit (or all hits)'''
    formatdb="formatdb"
    blastall="blastall"
    if os.path.isfile(os.path.join(ncbidir,"bin","formatdb")):
        formatdb     =os.path.join(ncbidir,"bin","formatdb")
        blastall     =os.path.join(ncbidir,"bin","blastall")
    elif os.path.isfile(os.path.join(ncbidir,"formatdb")):
        formatdb       =os.path.join(ncbidir,"formatdb")
        blastall       =os.path.join(ncbidir,"blastall")

    tmpdir="/tmp/"+os.getenv("USER")+'/'+os.path.basename(fastaA)+'_'+ \
        str(random.randrange(100000000))+os.path.basename(fastaB)+'/'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    tmpA=tmpdir+"fastaA"
    tmpB=tmpdir+"fastaB"
    shutil.copy(fastaA, tmpA)
    shutil.copy(fastaB, tmpB)

    fasta_listA,header_listA=split_fasta(tmpA)
    fasta_listB,header_listB=split_fasta(tmpB)
    
    fp=open(tmpA,'w')
    fp.write(''.join(fasta_listA))
    fp.close()
    
    os.system(formatdb+" -p T -o T -l "+tmpdir+"formatdb.log -i "+tmpA)

    # perform mapping by blastp 
    query=tmpdir+"query.fasta"
    map_dict=dict()
    for fastaB in fasta_listB:
        fp=open(query,'w')
        fp.write(fastaB)
        fp.close()
        command=blastall+" -p blastp -d "+tmpA+" -i "+query
        stdout,stderr=subprocess.Popen(command,
            shell=True, stdout = subprocess.PIPE).communicate()
        if not hit_txt in stdout:
            continue # No hits found
        hit_list=hit_pattern.findall(stdout.split(hit_txt)[1].split('>')[0])
        entryB=fastaB.splitlines()[0][1:]
        sys.stderr.write(entryB+'\n')
        if mapfirst:
            hit_list=hit_list[:1]
        for hit in hit_list:
            if hit in map_dict:
                map_dict[hit].append(entryB)
            else:
                map_dict[hit]=[entryB]
    
    map_txt='\n'.join(['%15s%5d '%(hit,len(map_dict[hit]))+' '.join(map_dict[hit]
        ) for hit in header_listA if hit in map_dict])+'\n'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    return map_txt


def map_library_phmmer(fastaA,fastaB,phmmer_path='phmmer',mapfirst=True):
    '''create phmmer mapping between two fasta sequence database "fastaA" 
    and "fastaB". Return mapping as text.
    "phmmer_path" is the path to phmmer
    "mapfirst" whether only take the first blastp hit (or all hits)'''

    tmpdir="/tmp/"+os.getenv("USER")+'/'+os.path.basename(fastaA)+'_'+ \
        str(random.randrange(100000000))+os.path.basename(fastaB)+'/'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    tmpA=tmpdir+os.path.basename(fastaA)
    tmpB=tmpdir+os.path.basename(fastaB)
    shutil.copy(fastaA, tmpA)
    shutil.copy(fastaB, tmpB)
    tblout=tmpdir+"tblout" # phmmmer tabular output

    fasta_listA,header_listA=split_fasta(tmpA)

    command=phmmer_path+" --noali --tblout "+tblout+' '+tmpB+' '+tmpA
    stdout,stderr=subprocess.Popen(command,
            shell=True, stdout = subprocess.PIPE).communicate()
    fp=open(tblout,'rU')
    tblout_txt=fp.read()
    fp.close()
    tblout_list=[line for line in tblout_txt.splitlines() if line.strip() \
        and not line.startswith('#')]

    # perform mapping by phmmer
    map_dict=dict()
    query_list=[]
    for line in tblout_txt.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        line=line.split()
        target=line[0]
        query=line[2]
        #evalue=float(line[4])
        if mapfirst and (query in query_list):
            continue
        else:
            query_list.append(query)

        if not target in map_dict:
            map_dict[target]=[query]
        else:
            map_dict[target].append(query)
    
    map_txt='\n'.join(['%15s%5d '%(hit,len(map_dict[hit]))+' '.join(map_dict[hit]
        ) for hit in header_listA if hit in map_dict])+'\n'
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    return map_txt

if __name__=="__main__":
    ncbidir=''
    mapfirst=True
    phmmer_path=''

    argv=[]
    for arg in sys.argv[1:]:
        if not arg.startswith('-'):
            argv.append(arg)
            continue
        if arg.startswith("-ncbidir="):
            ncbidir=arg[len("-ncbidir="):]
        elif arg.startswith("-mapfirst="):
            mapfirst=(arg[len("-mapfirst=")].lower()=="true")
        elif arg.startswith("-phmmer_path="):
            phmmer_path=arg[len("-phmmer_path="):]
        else:
            print >>sys.stderr, "ERROR! Unknown argument "+arg
            exit()

    if len(argv)<2:
        sys.stderr.write(docstring)
        exit()
    if phmmer_path and ncbidir:
        sys.stderr.write("ERROR! You cannot set both -phmmer_path and -ncbidir\n")
        exit()

    if not phmmer_path:
        map_txt=map_library(argv[0],argv[1],os.path.abspath(ncbidir),mapfirst)
    else:
        map_txt=map_library_phmmer(argv[0],argv[1],phmmer_path,mapfirst)

    if len(argv)>=3:
        fp=open(argv[2],'w')
        fp.write(map_txt)
        fp.close()
    else:
        sys.stdout.write(map_txt)
