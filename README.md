# HAST
Partition stLFR reads based on trio-binning algorithm using paternally unique markers.

## INSTALL

```
git clone https://github.com/BGI-Qingdao/HAST.git
cd HAST
make
```

or try the binary release https://github.com/BGI-Qingdao/HAST/releases/download/v1.0.0/HAST_release_v1.0.0_with_jellyfish.tar.gz

## USAGE

```
./HAST.sh [OPTION]

Trio-phase filial stLFR reads based on paternal and maternal NGS reads.

Options  :
        --paternal    paternal NGS reads file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --maternal    maternal NGS reads file in FASTA/FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --filial      filial stLFR reads file in FASTQ format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --thread      thread number.
                      [ optional, default 8 threads. ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      ( note: real memory used may be greater than this. )
                      [ optional, default 20GB. ]
        --jellyfish   jellyfish path.
                      [ optional, default jellyfish. ]
        --mer         mer-size
                      [ optional, default 21. ]
        --m-lower     maternal kmer frequency table will ignore kmers with count < m-lower.
                      [ optional, default 9. ]
        --m-upper     maternal kmer frequency table will ignore kmers with count > m-upper.
                      [ optional, default 33. ]
        --p-lower     paternal kmer frequency table will ignore kmers with count < p-lower.
                      [ optional, default 9. ]
        --p-upper     paternal kmer frequency table will ignore kmers with count > p-upper.
                      [ optional, default 33. ]
        --adaptor_f   forward adaptor
                      [ optional, default CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA. ]
        --adaptor_r   reverse adaptor
                      [ optional, default TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG. ]
        --auto_bounds automatically calcuate lower and upper bounds based on kmer analysis.
                      [ optional, default not trigger; no parameter. ]
                      ( note : if auto_bounds is on, it will overwrite --*-lower and --*-upper  ]
        --help        print this usage message.

Examples :
    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.fastq

    ./HAST.sh --paternal father.1.fastq.gz --paternal father.2.fastq.gz --maternal moter.fastq --filial son.fastq

    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.r1.fastq --filial son.r2.fastq

    ./HAST.sh --paternal father.fastq --maternal mother.fastq \
                     --filial son.r1.fastq --memory 20 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=33 --m-lower=9 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux
```

## Input and Output details examples

### Input

* scenario 01 ： full pipeine

        Parent's WGS sequences :  paternal.ngs.fastq  & maternal.ngs.fastq ;
        Children stLFR sequences : split_reads.1.fq.gz & split_reads.2.fq.gz ; NOTICE : stLFR raw reads are not supported , please split barcode first by https://github.com/BGI-Qingdao/stLFR_barcode_split.git

* scenario 02 ： classify_only

        Parent's unshare kmers          :  paternal.kmer  & maternal.kmer
        Children stLFR sequences : split_reads.1.fq.gz & split_reads.2.fq.gz ; NOTICE : stLFR raw reads are not supported , please split barcode first by https://github.com/BGI-Qingdao/stLFR_barcode_split.git

### Output
    
    maternal.split_reads.1.fq
    maternal.split_reads.2.fq
   
    paternal.split_reads.1.fq
    paternal.split_reads.2.fq
  
    homozygous.split_reads.1.fq
    homozygous.split_reads.2.fq
 
## How to run haplotype denovo after HAST

1. merge maternal & homozygous reads

```
    mkdir maternal
    cat maternal.split_reads.1.fq homozygous.split_reads.1.fq |gzip - >maternal/split_reads.1.fq.gz
    cat maternal.split_reads.2.fq homozygous.split_reads.2.fq |gzip - >maternal/split_reads.2.fq.gz
```
2. merge paternal & homozygous reads
```
    mkdir paternal
    cat paternal.split_reads.1.fq homozygous.split_reads.1.fq |gzip - >paternal/split_reads.1.fq.gz
    cat paternal.split_reads.2.fq homozygous.split_reads.2.fq |gzip - >paternal/split_reads.2.fq.gz
```

3. run stlfr2supernova pipeline independently by https://github.com/BGI-Qingdao/stlfr2supernova_pipeline#start-stlfr-clean

4. refine supernova assembly result by https://github.com/BGI-Qingdao/fabulous2.0

# run HSAT for 10X Genomics Linked Reads 

see https://github.com/BGI-Qingdao/HAST410GX

# run HAST for TGS reads

see https://github.com/BGI-Qingdao/HAST4TGS

____________________________
Enjoy !
