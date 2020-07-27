# HAST
Partition stLFR reads based on trio-binning algorithm using parentally unique markers.

## INSTALL

```
git clone https://github.com/BGI-Qingdao/HAST.git
cd HAST
make
```

## USAGE

An example script in HAST.sh

```
Usage   :
        HAST.sh <options>

Options :
        -h/--help       print this usage and exit ;

        --paternal      paternal ngs reads in fasta or fastq format ;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes (" ") ;

        --maternal      maternal ngs reads in fasta or fastq format ;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes (" ") ;
        --read1         filial read1 of stLFR reads in fastq format ;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes (" ") ;
                         * file name should contain "r1" ;

        --read2         filial read2 of stLFR reads;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes (" ") ;
                         * file name should contain "r2" ;

        --supernova     installed supernova path;

        --thread        max threads ;

Example :
        ./HAST.sh --supernova /home/software/supernova \
                  --paternal pa.wgs.fasta \
                  --maternal ma.wgs.fasta \
                  --read1 child.r1.fastq.gz \
                  --read2 child.r2.fastq.gz \

```

* one-step execution is allowed, but step-by-step running is recommended !!!

* Each step folder contains individual README.

# run HSAT for 10X Genomics Linked Reads 

see https://github.com/BGI-Qingdao/HAST410XG

# run HAST for TGS long reads (PacBio or Nanopore)

see https://github.com/BGI-Qingdao/HAST4TGS

____________________________
Enjoy !
