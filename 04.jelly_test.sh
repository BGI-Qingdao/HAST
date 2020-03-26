JELLYFISH='/ldfssz1/ST_OCEAN/USER/xumengyang/software/anaconda2/bin/jellyfish'

# JELLYFISH count , count all kmer's count.
# -m 21 : 21-mer.
# -s 20G : Initial hash size 20G ( human size ).
# -t 8 : 8 threads.
# -C : Count both strand, canonical representation.
$JELLYFISH count -m 21 -s 20G -t 8 -C -o  maternal_mer_counts.jf chr19.maternal.fastq
$JELLYFISH count -m 21 -s 20G -t 8 -C -o  paternal_mer_counts.jf chr19.paternal.fastq
# JELLYFISH dump , dump kmer-count as fasta format.
# -L 9 -U 33 : need kmer count within [9,33] interval.
$JELLYFISH dump -L 9 -U 33 maternal_mer_counts.jf -o maternal.mer.filter.fa
$JELLYFISH dump -L 9 -U 33 paternal_mer_counts.jf -o paternal.mer.filter.fa
$JELLYFISH dump maternal_mer_counts.jf            -o maternal.mer.fa
$JELLYFISH dump paternal_mer_counts.jf            -o paternal.mer.fa
# mixed 2 copy of maternal kmers + 1 copy of paternal kmers
cat maternal.mer.fa maternal.mer.fa paternal.mer.fa >mixed.fa
# JELLYFISH count
$JELLYFISH count -m 21 -s 20G -t 8 -C -o mixed_mer_counts.js mixed.fa
# JELLYFISH dump , dump kmer-count as fasta format.
# any kmer count as 1 belong to paternal only
# any kmer count as 2 belong to maternal only
$JELLYFISH dump -U 1 mixed_mer_counts.js      >paternal.mer.unique.fa
$JELLYFISH dump -L 2 -U 2 mixed_mer_counts.js >maternal.mer.unique.fa
# merge filtered & unique
cat paternal.mer.unique.fa paternal.mer.filter.fa > paternal_mixed.mer.fa
cat maternal.mer.unique.fa maternal.mer.filter.fa > maternal_mixed.mer.fa
# extract final kmers
$JELLYFISH count -m 21 -s 20G -t 8 -C -o paternal_mixed_mer_counts.js paternal_mixed.mer.fa
$JELLYFISH count -m 21 -s 20G -t 8 -C -o maternal_mixed_mer_counts.js maternal_mixed.mer.fa
# JELLYFISH dump , dump kmer-count by 2 column format , only first column is needed.
$JELLYFISH dump -t -c -L 2 -U 2 paternal_mixed_mer_counts.js | awk '{print $1}' >paternal.unique.filter.kmer
$JELLYFISH dump -t -c -L 2 -U 2 maternal_mixed_mer_counts.js | awk '{print $1}' >maternal.unique.filter.kmer
