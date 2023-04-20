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
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes  " "  ;

        --maternal      maternal ngs reads in fasta or fastq format ;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes  " "  ;
        --read1         filial read1 of stLFR reads in fastq format ;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes  " "  ;
                         * file name should contain "r1" ;

        --read2         filial read2 of stLFR reads;
                         * accept input file in gzip format, but the name should be ended by ".gz"
                         * accept multiple files, but should be seperated by whitespace; and all the filenames are within one pair of double quotes  " "  ;
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

# run HAST for 10X Genomics Linked Reads 

see https://github.com/BGI-Qingdao/HAST410XG

# run HAST for TGS long reads (PacBio or Nanopore)

see https://github.com/BGI-Qingdao/HAST4TGS


## Citing HAST
If you use HAST in your work, please cite:
Mengyang Xu, Lidong Guo, Xiao Du, Lei Li, Brock A Peters, Li Deng, Ou Wang, Fang Chen, Jun Wang, Zhesheng Jiang, Jinglin Han, Ming Ni, Huanming Yang, Xun Xu, Xin Liu, Jie Huang, Guangyi Fan, Accurate haplotype-resolved assembly reveals the origin of structural variants for human trios, Bioinformatics, Volume 37, Issue 15, 1 August 2021, Pages 2095–2102, https://doi.org/10.1093/bioinformatics/btab068

____________________________
Enjoy !
