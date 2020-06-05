#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./build_unshared_kmers.sh [OPTION]" 
    echo ""
    echo "Build parental unshared-kmers based on paternal and maternal NGS reads by jellyfish."
    echo ""
    echo "Options  :"
    echo "        --paternal    paternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --maternal    maternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --thread      thread number."
    echo "                      [ optional, default 8 threads. ]"
    echo "        --memory      x (GB) of memory to initial hash table by jellyfish."
    echo "                      ( note: real memory used may be greater than this. )"
    echo "                      [ optional, default 20GB. ]"
    echo "        --mer         mer-size"
    echo "                      [ optional, default 21. ]"
    echo "        --m-lower     maternal kmer frequency table will ignore kmers with count < m-lower."
    echo "                      [ optional, default 9. ]"
    echo "        --m-upper     maternal kmer frequency table will ignore kmers with count > m-upper."
    echo "                      [ optional, default 33. ]"
    echo "        --p-lower     paternal kmer frequency table will ignore kmers with count < p-lower."
    echo "                      [ optional, default 9. ]"
    echo "        --p-upper     paternal kmer frequency table will ignore kmers with count > p-upper."
    echo "                      [ optional, default 33. ]"
    echo "        --auto_bounds automatically calcuate lower and upper bounds based on kmer analysis."
    echo "                      [ optional, default not trigger; no parameter. ]"
    echo "                      ( note : if auto_bounds is on, it will overwrite --*-lower and --*-upper  ]"
    echo "        --help        print this usage message."
}

###############################################################################
# basic variables
###############################################################################
MER=21
CPU=8
MEMORY=10
PLOWER=9
PUPPER=33
MLOWER=9
MUPPER=33
PATERNAL=""
MATERNAL=""
AUTO_BOUNDS=0
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
        "--memory")
            MEMORY=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--m-lower")
            MLOWER=$2
            shift
            ;;
        "--m-upper")
            MUPPER=$2
            shift
            ;;
        "--p-lower")
            PLOWER=$2
            shift
            ;;
        "--p-upper")
            PUPPER=$2
            shift
            ;;
        "--mer")
            MER=$2
            shift
            ;;
        "--auto_bounds")
            AUTO_BOUNDS=1
            ;;
        "--paternal")
            PATERNAL=$2" "$PATERNAL
            shift
            ;;
        "--maternal")
            MATERNAL=$2" "$MATERNAL
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
echo "HAST starting with : "
echo "    paternal input : $PATERNAL"
echo "    maternal input : $MATERNAL"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower(maternal): $MLOWER"
echo "    upper(maternal): $MUPPER"
echo "    lower(paternal): $PLOWER"
echo "    upper(paternal): $PUPPER"
echo "    auto_bounds    : $AUTO_BOUNDS"
echo "build_unshared_kmer.sh in dir  : $SPATH"

JELLY=$SPATH"/jellyfish-linux"
ANALYSIS=$SPATH"/analysis_kmercount.sh"
if [[ ! -e $ANALYSIS  && $AUTO_BOUNDS == 1 ]] ; then
    echo "ERROR : \"$ANALYSIS\"  is missing. please download it from github. exit..."
    exit 1
fi
if [[ ! -e $JELLY ]] ; then
    echo "ERROR : \"$JELLY\"  is missing. please download it from github. exit..."
    exit 1
fi

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL  \
    $MER -lt 11 || \
    $MLOWER -lt 1 || $MUPPER -gt 100000000 || \
    $PLOWER -lt 1 || $PUPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
    exit 1
fi
for x in $MATERNAL $PATERNAL
do
   if [[ ! -e $x ]] ; then 
       echo "ERROR : input file \"$x\" is not exist ! exit ..."
       exit 1
   fi
done
date

###############################################################################
# extract paternal.unique.filter.mer & maternal.unique.filter.mer
###############################################################################
# count NGS reads
echo "extract unique mers by jellyfish ..."
if [[ ! -e "step_01_done" ]] ; then
    gz=0
    for fname in $MATERNAL
    do
        if [[ ${fname: -3} == ".gz"  ]] ; then
            if [[ $gz == 0 || $gz == 2 ]] ; then
                gz=2
            else
                echo "ERROR : please don't mixed gz input with non-gz input."
                exit 1;
            fi
        else
            if [[ $gz == 1 || $gz == 0 ]] ; then
                gz=1
            else
                echo "ERROR : please don't mixed gz input with non-gz input;"
                exit 1;
            fi
        fi
    done
    if [[ $gz == 2 ]] ; then
        zcat $MATERNAL | $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  maternal_mer_counts.jf /dev/fd/0 || exit 1
    else
        $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  maternal_mer_counts.jf  $MATERNAL || exit 1
    fi
    date >>"step_01_done"
else
    echo "skip kmer count of maternal because step_01_done file already exist ..."
fi

if [[ ! -e "step_02_done" ]] ; then
    gz=0
    for fname in $PATERNAL
    do
        if [[ ${fname: -3} == ".gz"  ]] ; then
            if [[ $gz == 0 || $gz == 2 ]] ; then
                gz=2
            else
                echo "ERROR : please don't mixed gz input with non-gz input."
                exit 1;
            fi
        else
            if [[ $gz == 1 || $gz == 0 ]] ; then
                gz=1
            else
                echo "ERROR : please don't mixed gz input with non-gz input;"
                exit 1;
            fi
        fi
    done
    if [[ $gz == 2 ]] ; then
        zcat $PATERNAL | $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  paternal_mer_counts.jf /dev/fd/0  || exit 1
    else
        $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o  paternal_mer_counts.jf $PATERNAL || exit 1
    fi
    date >>"step_02_done"
else
    echo "skip kmer count of paternal because step_02_done file already exist ..."
fi

# dump all mers
if [[ ! -e "step_03_done" ]] ; then
    $JELLY dump maternal_mer_counts.jf            -o maternal.mer.fa                    || exit 1
    date >>"step_03_done"
else
    echo "skip dump fa of maternal because step_03_done file already exist ..."
fi

if [[ ! -e "step_04_done" ]] ; then
    $JELLY dump paternal_mer_counts.jf            -o paternal.mer.fa                    || exit 1
    date >>"step_04_done"
else
    echo "skip dump fa of paternal because step_04_done file already exist ..."
fi

if [[ $AUTO_BOUNDS == 1 ]] ; then
    if  [[ ! -e "step_04.1_done" ]] ; then
        sh $ANALYSIS  || exit 1
        date >>"step_04.1_done"
    else 
        echo "skip kmer bounds analysis because step_04.1_done file already exist ..."
    fi
    MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
    PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
fi
echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "
# dump filter mers
if [[ ! -e "step_05_done" ]] ; then
    $JELLY dump -L $MLOWER -U $MUPPER maternal_mer_counts.jf -o maternal.mer.filter.fa || exit 1
    date >>"step_05_done"
else
    echo "skip dump fa of maternal filter  because step_05_done file already exist ..."
fi
if [[ ! -e "step_06_done" ]] ; then
    $JELLY dump -L $PLOWER -U $PUPPER paternal_mer_counts.jf -o paternal.mer.filter.fa || exit 1
    date >>"step_06_done"
else
    echo "skip dump fa of paternal filter  because step_06_done file already exist ..."
fi
# rm temporary files
rm -f maternal_mer_counts.jf paternal_mer_counts.jf
if [[ ! -e "step_07_done" ]] ; then
    # mix 1 copy of paternal mers and 2 copy of maternal mers and count p/maternal mixed mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o mixed_mer_counts.js  maternal.mer.fa maternal.mer.fa paternal.mer.fa  || exit 1
    # count==1 refer to paternal unique mers
    $JELLY dump -U 1 mixed_mer_counts.js          >paternal.mer.unique.fa  || exit 1
    # count==2 refer to maternal unique mers 
    $JELLY dump -L 2 -U 2 mixed_mer_counts.js     >maternal.mer.unique.fa || exit 1
    # rm temporary files
    rm  mixed_mer_counts.js
    date >>"step_07_done"
else
    echo "skip extract *aternal.mer.unique.fa  because step_07_done file already exist ..."
fi

if [[ ! -e "step_08_done" ]] ; then
    # count unique and filer mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o paternal_mixed_mer_counts.js paternal.mer.unique.fa paternal.mer.filter.fa || exit 1
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o maternal_mixed_mer_counts.js maternal.mer.unique.fa maternal.mer.filter.fa || exit 1
    # extrat both unique and filter mers
    $JELLY dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.mer || exit 1
    $JELLY dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.mer || exit 1
    # rm temporary files
    rm paternal_mixed_mer_counts.js
    rm maternal_mixed_mer_counts.js
    date >>"step_08_done"
else
    echo "skip extract *aternal.unique.filter.mer  because step_08_done file already exist ..."
fi

echo "final paternal unique kmer is : "
wc -l paternal.unique.filter.mer
echo "final maternal unique kmer is : "
wc -l maternal.unique.filter.mer

echo "extract unique mers done..."
date
