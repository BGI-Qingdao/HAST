#!/bin/bash

##
# change below variable as you wish
##
MIN_READPAIR_IN_BARCODE=1
THREADS=30
MEMORY=800


if [[ $# != 5 ]] ; then
    echo """
Usage   :
            ./assemble_by_supernova.sh supernova-path parental.r1.fq parental.r2.fq homo.r1.fq homo.r2.fq
"""
    exit 1;
fi

# step 00 : parse parameters and santity check
SUPERNOVA_PATH=$1
PS_1=$2
PS_2=$3
HO_1=$4
HO_2=$5
SUPERNOVA_WL=$SUPERNOVA_PATH/supernova-cs/*/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt
SUPERNOVA=$SUPERNOVA_PATH/supernova
SCRIPT_PATH=`dirname $0`

if [[ ! -e SUPERNOVA_WL ]] ; then 
    echo "supernova in not valid in $SUPERNOVA_PATH !!! ";
    echo "exit ..."
    exit 1;
fi

if [[ ! -e SUPERNOVA ]] ; then 
    echo "supernova in not valid in $SUPERNOVA_PATH !!! ";
    echo "exit ..."
    exit 1;
fi

if [[ ! -e $PS_1 || ! -e $PS_2 || ! -e $HO_1 || ! -e $HO_2 ]] ; then
    echo "input file invalid ... exit !!!"
    exit 1;
fi

# step 01 : merge parental.specific.fastq and homozygote.fastq

cat $PS_1 $HO_1 | gzip - >split_reads.1.fq.gz
cat $PS_2 $HO_2 | gzip - >split_reads.2.fq.gz

# step 02 : generate barcode_freq.txt
gzip -dc split_reads.1.fq.gz | awk '!(NR%4-1)' | awk -F '[# |]' '{print$2}' | awk -F '/' '{print $1}' | sort | uniq -c | awk '{printf("%s\t%s\n",$2,$1);}' > barcode_freq.txt

# mapping stLFR barcode into 10X Genomics
$SCRIPT_PATH/merge_barcodes.pl barcode_freq.txt  $SUPERNOVA_WL merge.txt $MIN_READPAIR_IN_BARCODE  1> merge_barcode.log  2>merge_barcode.err || exit 1

# format transfer from stLFR reass into 10X Genomics raw reads
$SCRIPT_PATH/fake_10x.pl  split_reads.1.fq.gz split_reads.2.fq.gz merge.txt >fake_10X.log 2>fake_10X.err || exit 1

# assembly by supernova
$SUPERNOVA run --id=haplotype --maxreads='all' --accept-extreme-coverage --fastqs=./ --localcores=$THREADS --localmem=$MEMORY --nopreflight >supernova_run.log 2>supernova_run.err || exit 1

$SUPERNOVA mkoutput --style=pseudohap2 --index --headers=full \
                --minsize=200 --asmdir=haplotype/outs/assembly/ \
                --outprefix=output > output.pshap2.log  2>&1
