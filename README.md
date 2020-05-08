# HAST
Partition stLFR reads based on trio-binning algorithm using paternally unique markers.

## INSTALL

```
git clone https://github.com/BGI-Qingdao/HAST.git
cd HAST
make
```

## USAGE

```
Usage    :
    ./HAST.sh [OPTION]

Trio-phase filial stLFR reads based on paternal and maternal NGS reads.

Options  :
        --paternal    paternal NGS reads file in FASTQ format.
                      ( note : gzip format is NOT supported. )
        --maternal    maternal NGS reads file in FASTQ format.
                      ( note : gzip format is NOT supported. )
        --filial      filial TGS reads file in FASTA format.
                      file in gzip format can be accepted, but filename must end by ".gz".
        --thread      thread num.
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
        --auto_bounds automatically calcuate lower and upper bounds based on kmer analysis.
                      [ optional, default not trigger; no parameter. ]
                      ( note : if auto_bounds is on, it will overwrite --*-lower and --*-upper  ]
        --help        print this usage message.

Examples :
    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.fastq

    ./HAST.sh --paternal father.fastq --maternal mother.fastq --filial son.r1.fastq --filial son.r2.fastq

    ./HAST.sh --paternal father.fastq --maternal mother.fastq \
                     --filial son.r1.fastq --memory 50 --thread 20 \
                     --mer 21 --p-lower=9 --p-upper=33 --m-lower=9 --p-upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux
```

Enjoy !
