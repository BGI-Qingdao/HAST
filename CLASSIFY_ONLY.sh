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
FILTER_FQ_BY_BARCODES_AWK=$SPATH"/filter_fq_by_barcodes.awk"

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
if [[ ! -e $FILTER_FQ_BY_BARCODES_AWK ]] ; then
    echo "ERROR : \"$FILTER_FQ_BY_BARCODES_AWK\"  is missing. please download it from github. exit..."
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
echo "extract unique barcode by classify ..."
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done
$CLASSIFY --hap0 $PATERNAL --hap1 $MATERNAL \
    --thread $CPU --weight1 1.04 $READ   >phased.barcodes 2>phased.log

awk '{if($2 == 0) print $1;}' phased.barcodes >paternal.unique.barcodes
echo "final paternal barcode :"
wc -l paternal.unique.barcodes
awk '{if($2 == 1) print $1;}' phased.barcodes >maternal.unique.barcodes
echo "final maternal barcodes"
wc -l maternal.unique.barcodes
awk '{if($2 == "-1") print $1;}' phased.barcodes >homozygous.unique.barcodes
echo "extract unique barcode done"
echo "final homozygous barcodes"
wc -l homozygous.unique.barcodes
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
date
echo "phase reads ..."
for x in $FILIAL
do
    name=`basename $x`
    if [[ ${name: -3} == ".gz" ]] ; then
        name=${name%%.gz}
        gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK maternal.unique.barcodes - >"maternal."$name
        gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK paternal.unique.barcodes - >"paternal."$name
        gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK homozygous.unique.barcodes - >"homozygous."$name
    else 
        awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  maternal.unique.barcodes $x >"maternal."$name
        awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  paternal.unique.barcodes $x >"paternal."$name
        awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  homozygous.unique.barcodes $x >"homozygous."$name
    fi
done
echo "phase reads done"
date
echo "__END__"
