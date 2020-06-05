

```
Usage    :
    ./classify_stlfr_reads.sh [OPTION]

Trio phase TGS reads use exist kmer datasets.

Options  :
        --paternal_mer paternal unique kmers
        --maternal_mer maternal unique kmers
        --filial       filial TGS reads file.
                       file in gzip format can be accepted, but filename must end by .gz.
        --thread       thread num.
                       [ optional, default 8 threads. ]
        --adaptor_f   forward adaptor
                      [ optional, default CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGA. ]
        --adaptor_r   reverse adaptor
                      [ optional, default TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGG. ]
        --help         print this usage message.

Examples :
    ./classify_stlfr_reads.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq

    ./classify_stlfr_reads.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --filial son2.fastq

    ./classify_stlfr_reads.sh --paternal_mer father.mers --maternal_mer mater.mers --filial son.fastq --thread 20

```
