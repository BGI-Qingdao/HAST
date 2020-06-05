#!/bin/bash

function usage(){
echo """
Usage   :
            ./assemble_by_supernova.sh <options>

Options :
                -h/--help           print this usage and exit;
                --supernova         supernova install path;
                --read1             read1 fastq;
                --read2             read2 fastq;
                --cpu               max threads to used by supernova;
                --memory            max memory(G) to used by supernova;
                --min_rp            minimax accept read-pair per barcode; barcode with read pair less-than min_rp will be ignored.

Example :
            ./assemble_by_supernova.sh --supernova /home/software/supernova \\
                --read1 paternal.r1.fastq  --read2 paternal.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \\
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \\
                --read2 paternal.r2.fastq --read2 homo.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \\
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \\
                --read2 paternal.r2.fastq --read2 homo.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \\
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \\
                --read2 paternal.r2.fastq --read2 homo.r2.fastq \\
                --cpu 30 --memory 500
"""
}

##
# change below variable as you wish
##
MIN_READPAIR_IN_BARCODE=1
THREADS=30
MEMORY=800
READ1=""
READ2=""
###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 ]] ; then 
    usage
    exit 0
fi
echo "CMD :$0 $*"
while [[ $# > 0 ]] 
do
    case $1 in
        "-h")
            usage
            exit 0
            ;;
        "--help")
            usage
            exit 0
            ;;
        "--supernova")
            SUPERNOVA_PATH=$2
            shift
            ;;
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            THREADS=$2
            shift
            ;;
        "--read2")
            READ2=$2" "$READ2
            shift 
            ;;
        "--read1")
            READ1=$2" "$READ1
            shift 
            ;;
        "--min_rp")
            MIN_READPAIR_IN_BARCODE=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

SUPERNOVA_WL=$SUPERNOVA_PATH/supernova-cs/*/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt
SUPERNOVA=$SUPERNOVA_PATH/supernova
SCRIPT_PATH=`dirname $0`

echo "INFO  : supernova in    : $SUPERNOVA_PATH"
echo "INFO  : all read1 files : $READ1"
echo "INFO  : all read2 files : $READ2"
echo "INFO  : max threads     : $CPU"
echo "INFO  : max memory      : $MEMORY"

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

for x in $READ1 $READ2
do 
    if [[ ! -e $x ]] ; then 
        echo "input file $x invalid ... exit !!!"
        exit 1;
    fi
done

# step 01 : merge parental.specific.fastq and homozygote.fastq

cat $READ1 | gzip  >split_reads.1.fq.gz
cat $READ2 | gzip  >split_reads.2.fq.gz

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
