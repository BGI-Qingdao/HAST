#!/bin/bash

function usage() {
    echo """
Usage   :
            mkoutput_by_fabulous2.0.sh <options>

Options :
            -h/--help           print this usage and exit;
            --assembly_path     supernova assembly path ; 
            --prefix            prefix of input & output files ;
                                please keep the prefix with assembly_by_supernova.sh ;
            --paternal_mer      paternal kmers ;
            --maternal_mer      maternal kmers ;
            --thread            max threads ;

Example :

./mkoutput_by_fabulous2.0.sh --assembly_path supernova_asmfolder \\
                                           --prefix  output \\
                                           --paternal_mer paternal.unique.filter.mer \\
                                           --maternal_mer maternal.unique.filter.mer
"""
}


# WARN      :
#               this is only a temporary script as example.
#               change below variables as you wish.
CPU=30
SUPERNOVA_OUTPUT_DIR=""
PATERNAL_KMERS=""
MATERNAL_KMERS=""
PREFIX="output"

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
        "--assembly_path")
            SUPERNOVA_OUTPUT_DIR=$2
            shift
            ;;
        "--prefix")
            PREFIX=$2
            shift
            ;;
        "--paternal_mer")
            PATERNAL_KMERS=$2
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--maternal_mer")
            MATERNAL_KMERS=$2
            shift
            ;;
        *)
            echo "invalid params : \"$1\" . exit ... "
            exit
        ;;
    esac
    shift
done

SCRIPT_PATH=`dirname $0`
WORKING_PATH=`pwd`
echo "LOG  : use scripts from $SCRIPT_PATH"
echo "LOG  : working dir $WORKING_PATH"
echo "LOG  : input supernova result in $SUPERNOVA_OUTPUT_DIR"
echo "LOG  : paternal kmers file is $PATERNAL_KMERS"
echo "LOG  : maternal kmers file is $MATERNAL_KMERS"
echo "LOG  : prefix is $PREFIX"
date

if [[ ! -e $SUPERNOVA_OUTPUT_DIR/$PREFIX.1.fasta || \
      ! -e $SUPERNOVA_OUTPUT_DIR/$PREFIX.2.fasta || \
      ! -e $SUPERNOVA_OUTPUT_DIR/$PREFIX.1.idx || \
      ! -e $SUPERNOVA_OUTPUT_DIR/$PREFIX.2.idx ]] ; then 
    echo "$SUPERNOVA_OUTPUT_DIR is not a valid path . exit ..."
    exit 1 ;
fi

# split hap1,hap2 from bubbles
$SCRIPT_PATH/bin/Split --fa_1 $SUPERNOVA_OUTPUT_DIR/${PREFIX}.1.fasta \
                       --fa_2 $SUPERNOVA_OUTPUT_DIR/${PREFIX}.2.fasta \
                       --idx_1 $SUPERNOVA_OUTPUT_DIR/${PREFIX}.1.idx \
                       --idx_2 $SUPERNOVA_OUTPUT_DIR/${PREFIX}.2.idx  \
                       --prefix ${PREFIX}

# haplotype bubbles based on paternal unique kmers
cat ${PREFIX}.phb.1.fa ${PREFIX}.phb.2.fa >  ${PREFIX}.phb.12.fa

$SCRIPT_PATH/bin/classify --hap $PATERNAL_KMERS  --hap $MATERNAL_KMERS \
                          --thread $CPU  --read ${PREFIX}.phb.12.fa \
                          --format fasta >phasing.out  2>phasing.log
cat phasing.out |grep Read |grep haplotype0 |awk '{print $2}' > ${PREFIX}.phb.12.father.idx
cat phasing.out |grep Read |grep haplotype1 |awk '{print $2}' > ${PREFIX}.phb.12.mother.idx
cat phasing.out |grep Read |grep ambiguous |awk '{print $2}' > ${PREFIX}.phb.12.ambiguous.idx

# cluster parental groups and output fasta
$SCRIPT_PATH/bin/MergePhaseResult  --prefix ${PREFIX}  \
              --father_ids ${PREFIX}.phb.12.father.idx \
              --mother_ids ${PREFIX}.phb.12.mother.idx \
              --homo_ids ${PREFIX}.phb.12.ambiguous.idx

# generate final fasta sequence
$SCRIPT_PATH/bin/GenSq --prefix ${PREFIX} 

echo "LOG  : ALL DONE"
echo "LOG  : final output is ${PREFIX}.father.fa and ${PREFIX}.mather.fa"
