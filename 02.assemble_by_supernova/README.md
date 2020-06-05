
```

Usage   :
            ./assemble_by_supernova.sh <options>

Options :
                -h/--help           print this usage and exit;
                --supernova         supernova install path;
                --read1             read1 fastq;
                --read2             read2 fastq;
                --cpu               max threads to used by supernova;
                --memory            max memory(G) to used by supernova;
                --min_rp            minimax accept read-pair per barcode; barcode with read pair less-than min_rp will be ignored;
                --prefix            prefix of output files;

Example :
            ./assemble_by_supernova.sh --supernova /home/software/supernova \
                --read1 paternal.r1.fastq  --read2 paternal.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \
                --read2 paternal.r2.fastq --read2 homo.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \
                --read2 paternal.r2.fastq --read2 homo.r2.fastq

            ./assemble_by_supernova.sh --supernova /home/software/supernova \
                --read1 paternal.r1.fastq --read1 homo.r1.fastq \
                --read2 paternal.r2.fastq --read2 homo.r2.fastq \
                --cpu 30 --memory 500

```
