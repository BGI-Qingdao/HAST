#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
echo """
Usage    :
    ./CLASSIFY_ONLY.sh [OPTION]

Trio phase TGS reads use exist kmer datasets.

Options  :
        --paternal_mer paternal unique kmers
        --maternal_mer maternal unique kmers
        --filial       filial TGS reads file.
                       file in gzip format can be accepted, but filename must end by ".gz".
        --thread       thread num.
                       [ optional, default 8 threads. ]
        --adaptor_f   forward adaptor
                      [ optional, default CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA. ]
        --adaptor_r   reverse adaptor
                      [ optional, default TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG. ]
        --help         print this usage message.

Examples :
    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --filial son2.fastq

    ./CLASSIFY_ONLY.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --thread 20
"""
}

###############################################################################
# basic variables 
###############################################################################
CPU=8
PATERNAL=""
MATERNAL=""
FILIAL=""
FORMAT='fastq'
SPATH=`dirname $0`
AF='CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA'
AR='TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG'

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
        "--adaptor_r")
            A_R=$2
            shift
            ;;
        "--adaptor_f")
            A_F=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--paternal_mer")
            PATERNAL=$2
            shift
            ;;
        "--maternal_mer")
            MATERNAL=$2
            shift
            ;;
        "--filial")
            FILIAL=$2" "$FILIAL
            shift 
            ;;
        "--format")
            FORMAT=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done
# print arguments
echo "CLASSIFY_ONLY starting with : "
echo "    paternal kmers : $PATERNAL"
echo "    maternal kmers : $MATERNAL"
echo "    filial input   : $FILIAL"
echo "    filial format  : $FORMAT"
echo "    thread         : $CPU "
echo "CLASSIFY_ONLY.sh in dir  : $SPATH"

CLASSIFY=$SPATH"/classify"
QUARTERING_FASTQ=$SPATH"/quartering_fastq.awk"

# sanity check
if [[ $CPU -lt 1 || -z $PATERNAL || -z $MATERNAL || -z $FILIAL ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
if [[ $FORMAT != 'fasta' && $FORMAT != 'fastq' ]] ; then 
    echo "ERROR : format invalid ... exit!!!"
    exit 1
fi
if [[ ! -e $CLASSIFY ]] ; then 
    echo "ERROR : please run \"make\" command in $SPATH before using this script! exit..."
    exit 1
fi
if [[ ! -e $QUARTERING_FASTQ ]] ; then
    echo "ERROR : \"$QUARTERING_FASTQ\"  is missing. please download it from github. exit..."
    exit 1
fi
for x in $MATERNAL $PATERNAL $FILIAL
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done
date
echo "__START__"
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done
if [[ ! -e "step_9_done" ]] ; then
    echo "extract unique barcode by classify ..."
    $CLASSIFY --hap0 $PATERNAL --hap1 $MATERNAL \
        --thread $CPU --weight1 1.04 $READ --adaptor_f $AF --adaptor_r $AR   >phased.barcodes 2>phased.log || exit 1
    date >>"step_9_done"
else
    echo "skip classify because step_9_done file already exist ..."
fi

if [[ ! -e "step_10_done" ]] ; then
    awk '{if($2 == 0) print $1;}' phased.barcodes >paternal.unique.barcodes || exit 1
    echo "final paternal barcode :"
    wc -l paternal.unique.barcodes
    awk '{if($2 == 1) print $1;}' phased.barcodes >maternal.unique.barcodes || exit 1
    echo "final maternal barcodes"
    wc -l maternal.unique.barcodes
    awk '{if($2 == "-1") print $1;}' phased.barcodes >homozygous.unique.barcodes || exit 1
    echo "extract unique barcode done"
    echo "final homozygous barcodes"
    wc -l homozygous.unique.barcodes
    date >>"step_10_done"
else
    echo "skip extract barcode because step_10_done file already exist ..."
fi
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
date
echo "phase reads ..."
if [[ ! -e "step_11_done" ]] ; then
    for x in $FILIAL
    do
        name=`basename $x`
        if [[ ${name: -3} == ".gz" ]] ; then
            name=${name%%.gz}
            gzip -dc $x | awk -v  prefix=$name -F '#|/' -f $QUARTERING_FASTQ paternal.unique.barcodes maternal.unique.barcodes homozygous.unique.barcodes  - || exit 1
        else
            awk -v  prefix=$name -F '#|/' -f $QUARTERING_FASTQ paternal.unique.barcodes maternal.unique.barcodes homozygous.unique.barcodes $x || exit 1
        fi
    done
    date >>"step_11_done"
else
    echo "skip extract reads because step_11_done file already exist ..."
    echo "basically , this means nothing changed by this running ! "
fi
wait
echo "phase reads done"
date
echo "__END__"
