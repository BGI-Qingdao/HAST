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

                                NOTE : the order of paternal_mer and maternal_mer is meaningful.
                                we will set the first [x]aterna_mer as primary output and second 
                                [x]aterna_mer as secondary output.

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
PREFER=''
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
            if [[ $PREFER == '' ]] ; then 
                PREFER='paternal'
            fi
            shift
            ;;
        "--thread")
            CPU=$2
            shift
            ;;
        "--maternal_mer")
            MATERNAL_KMERS=$2
            if [[ $PREFER == '' ]] ; then 
                PREFER='maternal'
            fi
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
echo "LOG  : prefer type is $PREFER"
echo "NOTE : you can swtich the order of --paternel_mer and --maternel_mer to swtich prefer type"
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
cat phasing.out |grep haplotype0 | awk '{printf("%s\t%s\n",$1,$3);}' > ${PREFIX}.phb.12.father.idx
cat phasing.out |grep haplotype1 | awk '{printf("%s\t%s\n",$1,$3);}'  > ${PREFIX}.phb.12.mother.idx
cat phasing.out |grep ambiguous  | awk '{printf("%s\t%s\n",$1,$3);}' > ${PREFIX}.phb.12.ambiguous.idx

# cluster parental groups and output fasta
$SCRIPT_PATH/bin/MergePhaseResult  --prefix ${PREFIX}  \
              --father_ids ${PREFIX}.phb.12.father.idx \
              --mother_ids ${PREFIX}.phb.12.mother.idx \
              --homo_ids ${PREFIX}.phb.12.ambiguous.idx

# generate final fasta sequence

if [[ $PREFER == "paternal" ]] ; then
    $SCRIPT_PATH/bin/GenSq --prefix ${PREFIX} --prefer 'pat'
else
    $SCRIPT_PATH/bin/GenSq --prefix ${PREFIX} --prefer 'mat'
fi

if [[ $PREFER == "paternal" ]] ; then
    ln -s ${PREFIX}.father.fa ${PREFIX}.primary.fa
    if [[ -e ${PREFIX}.mother.fa ]] ; then
        ln -s ${PREFIX}.mother.fa ${PREFIX}.secondary.fa
    fi
else 
    ln -s ${PREFIX}.mother.fa ${PREFIX}.primary.fa
    if [[ -e ${PREFIX}.father.fa ]] ; then
        ln -s ${PREFIX}.father.fa ${PREFIX}.secondary.fa
    fi
fi
echo "LOG  : ALL DONE"
echo "LOG  : final output is ${PREFIX}.primary.fa"
