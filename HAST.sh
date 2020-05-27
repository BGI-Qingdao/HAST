#!/bin/bash

###############################################################################
# usage function
###############################################################################
function usage(){
    echo "Usage    :"
    echo "    ./HAST.sh [OPTION]" 
    echo ""
    echo "Trio-phase filial stLFR reads based on paternal and maternal NGS reads."
    echo ""
    echo "Options  :"
    echo "        --paternal    paternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --maternal    maternal NGS reads file in FASTA/FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --filial      filial stLFR reads file in FASTQ format."
    echo "                      file in gzip format can be accepted, but filename must end by \".gz\"."
    echo "        --thread      thread number."
    echo "                      [ optional, default 8 threads. ]"
    echo "        --memory      x (GB) of memory to initial hash table by jellyfish."
    echo "                      ( note: real memory used may be greater than this. )"
    echo "                      [ optional, default 20GB. ]"
    echo "        --jellyfish   jellyfish path."
    echo "                      [ optional, default jellyfish. ]"
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
    echo "        --adaptor_f   forward adaptor"
    echo "                      [ optional, default CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA. ]"
    echo "        --adaptor_r   reverse adaptor"
    echo "                      [ optional, default TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG. ]"
    echo "        --auto_bounds automatically calcuate lower and upper bounds based on kmer analysis."
    echo "                      [ optional, default not trigger; no parameter. ]"
    echo "                      ( note : if auto_bounds is on, it will overwrite --*-lower and --*-upper  ]"
    echo "        --help        print this usage message."
    echo "        "
    echo "Examples :"
    echo "    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.fastq"
    echo ""
    echo "    ./HAST.sh --paternal father.1.fastq.gz --paternal father.2.fastq.gz --maternal moter.fastq --filial son.fastq"
    echo ""
    echo "    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.r1.fastq --filial son.r2.fastq"
    echo ""
    echo "    ./HAST.sh --paternal father.fastq --maternal mother.fastq \\"
    echo "                     --filial son.r1.fastq --memory 20 --thread 20 \\"
    echo "                     --mer 21 --p-lower=9 --p-upper=33 --m-lower=9 --p-upper=33 \\"
    echo "                     --jellyfish /home/software/jellyfish/jellyfish-linux"
}

###############################################################################
# basic variables 
###############################################################################
MER=21
JELLY=jellyfish
CPU=8
MEMORY=10
PLOWER=9
PUPPER=33
MLOWER=9
MUPPER=33
PATERNAL=""
MATERNAL=""
FILIAL=""
AUTO_BOUNDS=0
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
        "--filial")
            FILIAL=$2" "$FILIAL
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
echo "    filial input   : $FILIAL"
echo "    jellyfish      : $JELLY"
echo "    memory         : $MEMORY GB"
echo "    thread         : $CPU "
echo "    mer            : $MER "
echo "    lower(maternal): $MLOWER"
echo "    upper(maternal): $MUPPER"
echo "    lower(paternal): $PLOWER"
echo "    upper(paternal): $PUPPER"
echo "    auto_bounds    : $AUTO_BOUNDS"
echo "HAST.sh in dir  : $SPATH"

CLASSIFY=$SPATH"/classify"
FILTER_FQ_BY_BARCODES_AWK=$SPATH"/filter_fq_by_barcodes.awk"
ANALYSIS=$SPATH"/analysis_kmercount.sh"

# sanity check
if [[ $MEMORY -lt 1  || $CPU -lt 1 || \
    -z $PATERNAL || -z $MATERNAL || -z $FILIAL || \
    -z $JELLY  || $MER -lt 11 || \
    $MLOWER -lt 1 || $MUPPER -gt 100000000 || \
    $PLOWER -lt 1 || $PUPPER -gt 100000000 ]] ; then
    echo "ERROR : arguments invalid ... exit!!! "
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
rm maternal_mer_counts.jf paternal_mer_counts.jf
if [[ ! -e "step_07_done" ]] ; then
    # mix 1 copy of paternal mers and 2 copy of maternal mers
    cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa
    # count p/maternal mixed mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o mixed_mer_counts.js mixed.fa || exit1
    # count==1 refer to paternal unique mers
    $JELLY dump -U 1 mixed_mer_counts.js          >paternal.mer.unique.fa  || exit 1
    # count==2 refer to maternal unique mers 
    $JELLY dump -L 2 -U 2 mixed_mer_counts.js     >maternal.mer.unique.fa || exit 1
    # rm temporary files
    rm mixed.fa mixed_mer_counts.js
    date >>"step_07_done"
else
    echo "skip extract *aternal.mer.unique.fa  because step_07_done file already exist ..."
fi

if [[ ! -e "step_08_done" ]] ; then
    # mix unique mers and filter mers
    cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa
    cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa
    # count unique and filer mers
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa || exit 1
    $JELLY count -m $MER -s $MEMORY"G" -t $CPU -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa || exit 1
    # extrat both unique and filter mers
    $JELLY dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.mer || exit 1
    $JELLY dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.mer || exit 1
    # rm temporary files
    rm paternal_mixed.mer.fa paternal_mixed_mer_counts.js
    rm maternal_mixed.mer.fa maternal_mixed_mer_counts.js
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
###############################################################################
# phase filial barcode based on unique and filter mers of paternal and maternal
###############################################################################
echo "extract unique barcode by classify ..."
for x in $FILIAL
do 
    READ="$READ"" --read ""$x"
done

if [[ ! -e "step_09_done" ]] ; then
$CLASSIFY --hap0 paternal.unique.filter.mer --hap1 maternal.unique.filter.mer \
    --thread $CPU --weight1 1.04 $READ --adaptor_f $AF --adaptor_r $AR >phased.barcodes 2>phased.log || exit 1
    date >>"step_09_done"
else
    echo "skip run classify  because step_09_done file already exist ..."
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
            gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK maternal.unique.barcodes - >"maternal."$name &
            gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK paternal.unique.barcodes - >"paternal."$name &
            gzip -dc $x | awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK homozygous.unique.barcodes - >"homozygous."$name &
        else 
            awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  maternal.unique.barcodes $x >"maternal."$name & 
            awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  paternal.unique.barcodes $x >"paternal."$name & 
            awk  -F '#|/' -f $FILTER_FQ_BY_BARCODES_AWK  homozygous.unique.barcodes $x >"homozygous."$name &
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
