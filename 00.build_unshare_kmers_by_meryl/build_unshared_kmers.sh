#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./build_unshared_kmers.sh [OPTION]" 
    echo ""
    echo "Build parental unshared-kmers based on paternal and maternal NGS reads by meryl."
    echo ""
    echo "Options  :"
    echo "        --paternal    paternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --maternal    maternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --thread      thread number."
    echo "                      [ optional, default 8 threads. ]"
    echo "        --memory      x (GB) of memory to used by meryl."
    echo "                      [ optional, default 100GB. ]"
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
    echo "                      ( !!! WARN : default bounds is seted for 30X WGS reads , if your data is not close to 30X, please use your own bounds or simply open auto_bounds !!! ) "
    echo "        --help        print this usage message."
}

function get_real_path(){
    local vv=""
    for x in $1
    do
        real_p=`realpath $1`
        vv=$vv" "$real_p
    done
    echo $vv
}

###############################################################################
# basic variables
###############################################################################
MER=21
CPU=8
MEMORY=100
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
PATERNAL=`get_real_path $PATERNAL`
MATERNAL=`get_real_path $MATERNAL`

echo "script pat         : $SPATH"
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

MERYL_COUNT=$SPATH"/meryl/meryl.sh"
MERYL=$SPATH"/meryl/meryl"
ANALYSIS=$SPATH"/analysis_kmercount.sh"
if [[ ! -e $ANALYSIS  && $AUTO_BOUNDS == 1 ]] ; then
    echo "ERROR : \"$ANALYSIS\"  is missing. please download it from github. exit..."
    exit 1
fi
if [[ ! -e $MERYL ]] ; then
    echo "ERROR : \"$MERYL\"  is missing. please download it from github. exit..."
    exit 1
fi
if [[ ! -e $MERYL_COUNT ]] ; then
    echo "ERROR : \"$MERYL_COUNT\"  is missing. please download it from github. exit..."
    exit 1
fi

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 ||
    -z $PATERNAL || -z $MATERNAL  ||
    $MER -lt 11 ||
    $MLOWER -lt 1 || $MUPPER -gt 100000000 ||
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
echo "extract unique mers by meryl ..."
if [[ ! -e "step_01_done" ]] ; then
    mkdir -p maternal_meryl
    cd maternal_meryl
    $MERYL_COUNT $MATERNAL  || exit 1
    cd ..
    date >>"step_01_done"
else
    echo "skip kmer count of maternal because step_01_done file already exist ..."
fi

if [[ ! -e "step_02_done" ]] ; then
    mkdir -p paternal_meryl
    cd paternal_meryl
    $MERYL_COUNT $PATERNAL || exit 1
    cd ..
    date >>"step_02_done"
else
    echo "skip kmer count of paternal because step_02_done file already exist ..."
fi

if [[ $AUTO_BOUNDS == 1 ]] ; then
    if  [[ ! -e "step_02.1_done" ]] ; then
        sh $ANALYSIS  || exit 1
        date >>"step_02.1_done"
    else 
        echo "skip kmer bounds analysis because step_03.1_done file already exist ..."
    fi
    MLOWER=`grep LOWER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    MUPPER=`grep UPPER_INDEX maternal.bounds.txt| awk -F '=' '{print $2}'`
    PLOWER=`grep LOWER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
    PUPPER=`grep UPPER_INDEX paternal.bounds.txt| awk -F '=' '{print $2}'`
fi
echo "  the real used kmer-count bounds of maternal is [ $MLOWER , $MUPPER ] "
echo "  the real used kmer-count bounds of paternal is [ $PLOWER , $PUPPER ] "
# dump filter mers
if [[ ! -e "step_03_done" ]] ; then
    $MERYL difference maternal_meryl/reads.final.meryl paternal_meryl/reads.final.meryl output mat_not_pat.meryl  || exit 1
    date >>"step_03_done"
else
    echo "skip step 03"
fi
if [[ ! -e "step_04_done" ]] ; then
    $MERYL difference paternal_meryl/reads.final.meryl maternal_meryl/reads.final.meryl output pat_not_mat.meryl || exit 1
    date >>"step_04_done"
else
    echo "skip step 04"
fi
if [[ ! -e "step_05_done" ]] ; then
    MLOWER=$(($MLOWER-1))
    MUPPER=$(($MUPPER+1))
    $MERYL greater-than $MLOWER  mat_not_pat.meryl output 'mat_not_pat_gt'$MLOWER'.meryl' || exit  1
    $MERYL less-than $MUPPER  'mat_not_pat_gt'$MLOWER'.meryl' output 'mat_not_pat_gt'$MLOWER'_lt'$MUPPER'.meryl' || exit 1
    date >>"step_05_done"
else
    echo "skip step 05"
fi

if [[ ! -e "step_06_done" ]] ; then
    PLOWER=$(($PLOWER-1))
    PUPPER=$(($PUPPER+1))
    $MERYL greater-than $PLOWER  pat_not_mat.meryl output 'pat_not_mat_gt'$PLOWER'.meryl' || exit  1
    $MERYL less-than $PUPPER  'pat_not_mat_gt'$PLOWER'.meryl' output 'pat_not_mat_gt'$PLOWER'_lt'$PUPPER'.meryl' || exit 1
    date >>"step_06_done"
else
    echo "skip step 06"
fi
if [[ ! -e "step_07_done" ]] ; then
    $MERYL print  'mat_not_pat_gt'$MLOWER'_lt'$MUPPER'.meryl' | awk '{print $1}' > maternal.unique.filter.mer || exit 1
    $MERYL print  'pat_not_mat_gt'$PLOWER'_lt'$PUPPER'.meryl' | awk '{print $1}' > paternal.unique.filter.mer || exit 1
    date >>"step_07_done"
else
    echo "skip extract *aternal.unique.filter.mer  because step_07_done file already exist ..."
fi

echo "final paternal unique kmer is : "
wc -l paternal.unique.filter.mer
echo "final maternal unique kmer is : "
wc -l maternal.unique.filter.mer

echo "extract unique mers done..."
date
