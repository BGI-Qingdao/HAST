#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./trioSLR.sh [OPTION]" 
    echo ""
    echo "Trio-phase filial stLFR reads based on paternal NGS reads and maternal NGS reads."
    echo ""
    echo "Options  :"
    echo "        --paternal    paternal NGS reads file in fastq format."
    echo "        --maternal    maternal NGS reads file in fastq format."
    echo "        --filial      filial stLFR reads file in fastq format."
    echo "        --thread      threads num."
    echo "                      [ optioal , default 8 thread ]"
    echo "        --memory      x (GB) of memory to be used by jellyfish."
    echo "                      [ optioal , default 20GB ]"
    echo "        --jellyfish   jellyfish path."
    echo "                      [ optioal , default jellyfish ]"
    echo "        --mer         mer-size"
    echo "                      [ optioal , default 21 ]"
    echo "        --lower       ignore mer with cout < lower."
    echo "                      [ optioal , default 9 ]"
    echo "        --upper       ignore mer with cout > upper."
    echo "                      [ optioal , default 33 ]"
    echo "        --help        print this usage message."
    echo "        "
    echo "Examples :"
    echo "    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq"
    echo ""
    echo "    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq --filial son.r1.fastq --filial son.r2.fastq"
    echo ""
    echo "    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq \\"
    echo "                     --filial son.r1.fastq --memory 20 --thread 20 \\"
    echo "                     --mer 21 --lower=9 --upper=33 \\"
    echo "                     --jellyfish /home/software/jellyfish/jellyfish-linux"
}

###############################################################################
# basic variables 
###############################################################################
MER=21
JELLY=jellyfish
CPU=8
MEMORY=10
LOWER=9
UPPER=33
PATERNAL=""
MATERNAL=""
FILIAL=""

###############################################################################
# parse arguments
###############################################################################
if [[ $# == 0 ]] ; then 
    usage
    exit 0
fi
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
        "--jellyfish")
            JELLY=$2
            shift
            ;;
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--lower")
            LOWER=$2
            shift
            ;;
        "--upper")
            UPPER=$2
            shift
            ;;
        "--mer")
            MER=$2
            shift
            ;;
        "--paternal")
            PATERNAL=$2
            shift
            ;;
        "--maternal")
            MATERNAL=$2
            shift
            ;;
        "--filial")
            FILIAL=$2" "$FILIAL
    esac
    shift
done
# print arguments
echo "trioSLR starting with : "
echo "    paternal input : $PATERNAL"
echo "    maternal input : $MATERNAL"
echo "    filial input   : $FILIAL"
echo "    jellyfish      : $JELLY"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower          : $LOWER"
echo "    upper          : $UPPER"
# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL || -z $FILIAL || \
    -z $JELLY  || $MER -lt 11 || \
    $LOWER -lt 1 || $UPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
###############################################################################
# extract paternal.unique.filter.mer & maternal.unique.filter.mer
###############################################################################
# count NGS reads
$JELLYFISH count -m $MER -s $MEMORY"G" -t $CPU -C -o  maternal_mer_counts.jf $MATERNAL
$JELLYFISH count -m $MER -s $MEMORY"G" -t $CPU -C -o  paternal_mer_counts.jf $PATERNAL
# dump filter mers
$JELLYFISH dump -L $LOWER -U $UPPER maternal_mer_counts.jf -o maternal.mer.filter.fa
$JELLYFISH dump -L $LOWER -U $UPPER paternal_mer_counts.jf -o paternal.mer.filter.fa
# dump all mers
$JELLYFISH dump maternal_mer_counts.jf            -o maternal.mer.fa
$JELLYFISH dump paternal_mer_counts.jf            -o paternal.mer.fa
# mix 1 copy of paternal mers and 2 copy of maternal mers
cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa
# count p/maternal mixed mers
$JELLYFISH count -m $MER -s $MEMORY"G" -t $CPU -C -o mixed_mer_counts.js mixed.fa
# count==1 refer to paternal unique mers
$JELLYFISH dump -U 1 mixed_mer_counts.js          >paternal.mer.unique.fa
# count==2 refer to maternal unique mers
$JELLYFISH dump -L 2 -U 2 mixed_mer_counts.js     >maternal.mer.unique.fa
# mix unique mers and filter mers
cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa
cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa
# count unique and filer mers
$JELLYFISH count -m $MER -s $MEMORY"G" -t $CPU -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa
$JELLYFISH count -m $MER -s $MEMORY"G" -t $CPU -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa
# extrat both unique and filter mers
$JELLYFISH dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.mer
$JELLYFISH dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.mer

###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done
./classify --hap0 paternal.unique.filter.mer --hap1 maternal.unique.filter.mer \
    --thread $CPU $READ >phasing.barcode 2>phasing.log

awk '{if($2 == 0) print $1;}' phased.barcodes >paternal.unique.barcodes
awk '{if($2 == 1) print $1;}' phased.barcodes >maternal.unique.barcodes
awk '{if($2 == "-1") print $1;}' phased.barcodes >homozygous.unique.barcodes
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
for x in $FILIAL
do
    name=`basename $x`
    awk  -F '#|/' -f filter_fq_by_barcodes.awk  maternal.unique.barcodes $x >"maternal."$name
    awk  -F '#|/' -f filter_fq_by_barcodes.awk  paternal.unique.barcodes $x >"paternal."$name
    awk  -F '#|/' -f filter_fq_by_barcodes.awk  homozygous.unique.barcodes $x >"homozygous."$name
done
