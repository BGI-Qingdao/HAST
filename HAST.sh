#!/bin/bash

PATERNAL='xxx.paternal.fasta'
MATERNAL='xxx.maternal.fasta'
C1='child.r1.fastq'
C2='child.r2.fastq'
SUPERNOVA='xxx/supernova'

SCRIPT_PATH=`dirname $0`

STEP0=$SCRIPT_PATH/'00.build_unshare_kmers_by_jellyfish/build_unshared_kmers.sh'
STEP1=$SCRIPT_PATH/'01.classify_stlfr_reads/classify_stlfr_reads.sh'
STEP2=$SCRIPT_PATH/'02.assemble_by_supernova/assemble_by_supernova.sh'
STEP3=$SCRIPT_PATH/'mkoutput_by_fabulous2.0.sh'
# step 00 build unshare kmers from paternal ngs reads

mkdir '00.build_kmers'
cd '00.build_kmers' 
$STEP0 --paternal $PATERNAL --maternal $MATERNAL
cd ..
mv '00.build_kmers/paternal.unique.filter.mer' ./
mv '00.build_kmers/maternal.unique.filter.mer' ./

# step 01 : classify filial stLFR reads

mkdir '01.classify_reads'
cd '01.classify_reads' 
$STEP1 --paternal_mer ../paternal.unique.filter.mer \
                      --maternal_mer ../maternal.unique.filter.mer \
                      --filter $C1 --filter $C2 >log 2>err
cd ..
mv '01.classify_reads/*.paternal.fastq' ./
mv '01.classify_reads/*.maternal.fastq' ./
mv '00.build_kmers/*.homozygote.fastq' ./
mv '00.build_kmers/*.nobarcode.fastq' ./


# step 02 : assembly by supernova
mkdir '02.maternal_assembly'
cd '02.maternal_assembly'
$STEP2 --supernova $SUPERNOVA_PATH --read1 ../*1*.maternal.fastq \
                                   --read1 ../*1*.nobarcode.fastq \
                                   --read2 ../*2*.maternal.fastq \
                                   --read2 ../*2*.nobarcode.fastq \
                                   --prefix output \
                                   >log 2>err
cd ..

mkdir '02.paternal_assembly'
cd '02.paternal_assembly'
$STEP2 --supernova $SUPERNOVA_PATH --read1 ../*1*.paternal.fastq \
                                   --read1 ../*1*.nobarcode.fastq \
                                   --read2 ../*2*.paternal.fastq \
                                   --read2 ../*2*.nobarcode.fastq \
                                   --prefix output \
                                   >log 2>err
cd ..

# step 03 : assembly by supernova
mkdir '03.maternal_output'
cd '03.maternal_output'
$STEP3 --supernova_path '../02.maternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
cd ..

mkdir '03.paternal_output'
cd '03.paternal_output'
$STEP3 --supernova_path '../02.paternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
cd ..

echo 'All done'
echo 'Final result is : 03.maternal_output/output.mather.fa and  03.paternal_output/output.father.fa'
echo "Bye"
