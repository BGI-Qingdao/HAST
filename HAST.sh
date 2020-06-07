#!/bin/bash

function usage() {
echo """
Usage   :
        HAST.sh <options>

Options :
        -h/--help       print this usage and exit ;

        --paternal      paternal ngs reads in fasta or fastq format ;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by \" \";

        --maternal      maternal ngs reads in fasta or fastq format ;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by \" \";

        --read1         filial read1 of stLFR reads in fastq format ;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by \" \";
                         * file name must contain \"r1\" , otherwise error will occur ;

        --read2         filial read2 of stLFR reads;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by \" \";
                         * file name must contain \"r2\" , otherwise error will occur ;

        --supernova     supernova install path;

        --thread        max thread ;

Example :
        ./HAST.sh --supernova /home/software/supernova \\
                  --paternal pa.wgs.fasta \\
                  --maternal ma.wgs.fasta \\
                  --read1 child.r1.fastq.gz \\
                  --read2 child.r2.fastq.gz \\

"""
}

PATERNAL=''
MATERNAL=''
C1=''
C2=''
SUPERNOVA='xxx/supernova'
CPU=8
MEMORY=800

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
        "--supernova")
            SUPERNOVA=$2
            shift
            ;;
        "--read2")
            C2=$2" "$C2
            shift
            ;;
        "--read1")
            C1=$2" "$C1
            shift
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

SCRIPT_PATH=`dirname $0`
STEP0=$SCRIPT_PATH/'00.build_unshare_kmers_by_jellyfish/build_unshared_kmers.sh'
STEP1=$SCRIPT_PATH/'01.classify_stlfr_reads/classify_stlfr_reads.sh'
STEP2=$SCRIPT_PATH/'02.assemble_by_supernova/assemble_by_supernova.sh'
STEP3=$SCRIPT_PATH/'03.mkoutput_by_fabulous2.0/mkoutput_by_fabulous2.0.sh'
# step 00 build unshare kmers from paternal ngs reads

if [[ ! -e $STEP0  || 
      ! -e $STEP1  ||
      ! -e $STEP2  ||
      ! -e $STEP3 ]] ; then 
    echo "FATAL :   required script is missing. exit ... "
    exit 1
fi

if [[ ! -e $SUPERNOVA/supernova ]] ; then 
    echo "$SUPERNOVA is not valid supernova path. exit ... "
    exit 1
fi

mkdir '00.build_kmers'
cd '00.build_kmers'
# log command first
echo """
$STEP0 --paternal " $PATERNAL " --maternal " $MATERNAL " \
       --thread $CPU --memory $MEMORY 
"""
$STEP0 --paternal " $PATERNAL " --maternal " $MATERNAL " \
       --thread $CPU --memory $MEMORY 
cd ..
mv '00.build_kmers/paternal.unique.filter.mer' ./
mv '00.build_kmers/maternal.unique.filter.mer' ./

# step 01 : classify filial stLFR reads

mkdir '01.classify_reads'
cd '01.classify_reads' 
echo """
$STEP1 --paternal_mer ../paternal.unique.filter.mer \
                      --maternal_mer ../maternal.unique.filter.mer \
                      --auto_bounds \
                      --filter " $C1 " \
                      --filter " $C2 " >log 2>err
"""
$STEP1 --paternal_mer ../paternal.unique.filter.mer \
                      --maternal_mer ../maternal.unique.filter.mer \
                      --auto_bounds \
                      --filter " $C1 " \
                      --filter " $C2 " >log 2>err
cd ..
mv '01.classify_reads/*.paternal.fastq' ./
mv '01.classify_reads/*.maternal.fastq' ./
mv '00.build_kmers/*.homozygote.fastq' ./
mv '00.build_kmers/*.nobarcode.fastq' ./


# step 02 : assembly by supernova
mkdir '02.maternal_assembly'
cd '02.maternal_assembly'
echo """
$STEP2 --supernova $SUPERNOVA_PATH --read1 ../*r1*.maternal.fastq \
                                   --read1 ../*r1*.nobarcode.fastq \
                                   --read2 ../*r2*.maternal.fastq \
                                   --read2 ../*r2*.nobarcode.fastq \
                                   --prefix output \
                                   >log 2>err
"""
$STEP2 --supernova $SUPERNOVA_PATH --read1 ../*r1*.maternal.fastq \
                                   --read1 ../*r1*.nobarcode.fastq \
                                   --read2 ../*r2*.maternal.fastq \
                                   --read2 ../*r2*.nobarcode.fastq \
                                   --prefix output \
                                   >log 2>err
# above codes shows why read1 must contain r1 and read2 must contain r2
cd ..

mkdir '02.paternal_assembly'
cd '02.paternal_assembly'

echo """
$STEP2 --supernova $SUPERNOVA_PATH --read1 ../*1*.paternal.fastq \
                                   --read1 ../*1*.nobarcode.fastq \
                                   --read2 ../*2*.paternal.fastq \
                                   --read2 ../*2*.nobarcode.fastq \
                                   --prefix output \
                                   >log 2>err
"""
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
echo """
$STEP3 --supernova_path '../02.maternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
"""
$STEP3 --supernova_path '../02.maternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
cd ..

mkdir '03.paternal_output'
cd '03.paternal_output'
echo """
$STEP3 --supernova_path '../02.paternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
"""
$STEP3 --supernova_path '../02.paternal_assembly' \
       --paternal_mer ../paternal.unique.filter.mer \
       --maternal_mer ../maternal.unique.filter.mer \
       --prefix output > log 2>err
cd ..

echo 'All done'
echo 'Final result is : 03.maternal_output/output.mather.fa and  03.paternal_output/output.father.fa'
echo "Bye"
