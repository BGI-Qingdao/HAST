#!/bin/bash

JELLYFISH='/ldfssz1/ST_OCEAN/USER/xumengyang/software/anaconda2/bin/jellyfish'
CPU=10
MEMORY="20G"
KMER_SIZE=21
MIN_KMER_COUNT=9   # >=MIN_KMER_COUNT
MAX_KMER_COUNT=33  # <=MAX_KMER_COUNT
TEMP_PREFIX="temp"
PATERNEL_FQ='chr19.paternal.fastq'
MATERNEL_FQ='chr19.maternal.fastq'
CHILD_READ1_FQ="chr19.child.fastq"
CHILD_READ2_FQ=""
####################################################################
#  STEP 01 : find paternal.unique.kmer & maternal.unique.kmer
####################################################################
JELLYFISH='/ldfssz1/ST_OCEAN/USER/xumengyang/software/anaconda2/bin/jellyfish'

# JELLYFISH count , count all kmer's count.
# -m 21 : 21-mer.
# -s 20G : Initial hash size 20G ( human size ).
# -t 8 : 8 threads.
# -C : Count both strand, canonical representation.
$JELLYFISH count -m 21 -s 20G -t 8 -C -o  maternal_mer_counts.jf chr19.maternal.fastq
$JELLYFISH count -m 21 -s 20G -t 8 -C -o  paternal_mer_counts.jf chr19.paternal.fastq
# JELLYFISH dump , dump kmer-count as fasta format.
# -L 9 -U 33 : need kmer count within [9,33] interval.
$JELLYFISH dump -L 9 -U 33 maternal_mer_counts.jf -o maternal.mer.filter.fa
$JELLYFISH dump -L 9 -U 33 paternal_mer_counts.jf -o paternal.mer.filter.fa
$JELLYFISH dump maternal_mer_counts.jf            -o maternal.mer.fa
$JELLYFISH dump paternal_mer_counts.jf            -o paternal.mer.fa
# mixed 2 copy of maternal kmers + 1 copy of paternal kmers
cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa
# JELLYFISH count
$JELLYFISH count -m 21 -s 20G -t 8 -C -o mixed_mer_counts.js mixed.fa
# JELLYFISH dump , dump kmer-count as fasta format.
# any kmer count as 1 belong to paternal only
# any kmer count as 2 belong to maternal only
$JELLYFISH dump -U 1 mixed_mer_counts.js      >paternal.mer.unique.fa
$JELLYFISH dump -L 2 -U 2 mixed_mer_counts.js >maternal.mer.unique.fa
# merge filtered & unique
cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa
cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa
# extract final kmers
$JELLYFISH count -m 21 -s 20G -t 8 -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa
$JELLYFISH count -m 21 -s 20G -t 8 -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa
# JELLYFISH dump , dump kmer-count by 2 column format , only first column is needed.
$JELLYFISH dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.kmer
$JELLYFISH dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.kmer

####################################################################
#  STEP 02: classify all barcodes
####################################################################
if [[ $CHILD_READ2_FQ == "" ]] ;then
    ./classify --hap0 paternal.unique.filter.kmer --hap1 maternal.unique.filter.kmer \
        --read1 chr19.child.fa --thread $CPU >phased.barcodes
else 
    ./classify --hap0 paternal.unique.filter.kmer --hap1 maternal.unique.filter.kmer \
        --read1 chr19.child.fa \
        --read2 chr19.child.fa --thread $CPU >phased.barcodes
fi
awk '{if($2 == 0) print $1;}' phased.barcodes >paternal_unique.barcodes
awk '{if($2 == 1) print $1;}' phased.barcodes >maternal_unique.barcodes
awk '{if($2 == "-1") print $1;}' phased.barcodes >homozygous_unique.barcodes
####################################################################
#  STEP 03: classify all reads
####################################################################
awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal_unique.barcodes $CHILD_READ1 > maternal.1.fastq
awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal_unique.barcodes $CHILD_READ1 > paternal.1.fastq
awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' homozygous_unique.barcodes $CHILD_READ1 > homozygous.1.fastq

if [[ $CHILD_READ2 == "" ]] ;then
    awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' maternal_unique.barcodes $CHILD_READ2 > maternal.2.fastq
    awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' paternal_unique.barcodes $CHILD_READ2 > paternal.2.fastq
    awk  -F '#|/'  ' {if( FILENAME == ARGV[1] ) { s[$1]=1} else { if(FNR%4==1&&NF>1){ if ($2 in s ){ print $0 ; c=1;} else {c=0} } else { if(c==1) { print $0 ; }  } } }' homozygous_unique.barcodes $CHILD_READ2 > homozygous.2.fastq
fi
rm -rf $TEMP_PREFIX".*"