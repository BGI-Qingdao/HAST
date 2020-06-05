#!/bin/bash
#
# WARN      :
#               this is only a temporary script as example.
#               change below variables as you wish.
CPU=30

# step 0. check parameters or print usage
if [[ $# != 3 ]] ; then 
    echo """Usage     :
            ./example.sh  dir-to-supernova-outs paternal.kmers maternal.kmers
Example   :
            ./example.sh  /home/project/supernova_test/human paternal.kmers maternal.kmers
    """
    exit 1
fi

SCRIPT_PATH=`dirname $0`
WORKING_PATH=`pwd`
SUPERNOVA_OUTPUT_DIR=$1
PATERNAL_KMERS=$2
MATERNAL_KMERS=$3
echo "LOG  : use scripts from $SCRIPT_PATH"
echo "LOG  : working dir $WORKING_PATH"
echo "LOG  : input supernova result in $SUPERNOVA_OUTPUT_DIR"
echo "LOG  : paternal kmers file is $PATERNAL_KMERS"
echo "LOG  : maternal kmers file is $MATERNAL_KMERS"
echo "WARN : please guarantee supernova is usable"
date

PREFIX='output'
cp $SUPERNOVA_OUTPUT_DIR/output* ./

# output supernova pse2
gunzip ${PREFIX}.1.fasta.gz
gunzip ${PREFIX}.2.fasta.gz

# split hap1,hap2 from bubbles
$SCRIPT_PATH/bin/Split --fa_1 ${PREFIX}.1.fasta --fa_2 ${PREFIX}.2.fasta \
                       --idx_1 ${PREFIX}.1.idx  --idx_2 ${PREFIX}.2.idx  \
                       --parefix ${PREFIX}

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

$SCRIPT_PATH/bin/GenSq --prefix ${PREFIX} ${PREFIX}.homo.fa \
                      ${PREFIX}.phb.1.fa ${PREFIX}.phb.2.fa \
                      ${PREFIX}.merge.father.ids \
                      ${PREFIX}.merge.mother.ids

echo "LOG  : ALL DONE"
echo "LOG  : final output is ${PREFIX}.phb.1.fa and ${PREFIX}.phb.2.fa"
