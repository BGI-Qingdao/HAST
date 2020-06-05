# HAST
Partition stLFR reads based on trio-binning algorithm using paternally unique markers.

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
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by " ";

        --maternal      maternal ngs reads in fasta or fastq format ;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by " ";

        --read1         filial read1 of stLFR reads in fastq format ;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by " ";
                         * file name must contain "r1" , otherwise error will occur ;

        --read2         filial read2 of stLFR reads;
                         * if file is gzip format, it must be ended by .gz
                         * multi-file-name should be sperated by whitespace and include by " ";
                         * file name must contain "r2" , otherwise error will occur ;

        --supernova     supernova install path;

        --thread        max thread ;

Example :
        ./HAST.sh --supernova /home/software/supernova \
                  --paternal pa.wgs.fasta \
                  --maternal ma.wgs.fasta \
                  --read1 child.r1.fastq.gz \
                  --read2 child.r2.fastq.gz \

```

* However ,to get more flexible, I recommand you to run this step by step  !!!

* There is separeta README in each step folder.

# run HSAT for 10X Genomics Linked Reads 

see https://github.com/BGI-Qingdao/HAST410GX

# run HAST for TGS reads

see https://github.com/BGI-Qingdao/HAST4TGS

____________________________
Enjoy !
