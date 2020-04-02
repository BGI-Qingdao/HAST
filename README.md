# trioSLR
trio-phase stLFR reads based on paternal and maternal NGS reads.

## INSTALl

```
git clone https://github.com/BGI-Qingdao/trioSLR.git
cd trioSLR
make
```

## USAGE

```
Usage    :
    ./trioSLR.sh [OPTION]

Trio-phase filial stLFR reads based on paternal NGS reads and maternal NGS reads.

Options  :
        --paternal    paternal NGS reads file in fastq format.
        --maternal    maternal NGS reads file in fastq format.
        --filial      filial stLFR reads file in fastq format.
                      file in gzip format is accepted, but filename must end by .gz
        --thread      threads num.
                      [ optioal , default 8 thread ]
        --memory      x (GB) of memory to initial hash table by jellyfish.
                      (noted: real memory used maybe greater than this )
                      [ optioal , default 20GB ]
        --jellyfish   jellyfish path.
                      [ optioal , default jellyfish ]
        --mer         mer-size
                      [ optioal , default 21 ]
        --lower       ignore mer with cout < lower.
                      [ optioal , default 9 ]
        --upper       ignore mer with cout > upper.
                      [ optioal , default 33 ]
        --help        print this usage message.

Examples :
    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq --filial son.fastq

    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq --filial son.r1.fastq --filial son.r2.fastq

    ./trioSLR.sh --paternal father.fastq --maternal mater.fastq \
                     --filial son.r1.fastq --memory 20 --thread 20 \
                     --mer 21 --lower=9 --upper=33 \
                     --jellyfish /home/software/jellyfish/jellyfish-linux
```

Enjoy !
