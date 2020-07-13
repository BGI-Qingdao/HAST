#!/bin/bash

path=`dirname $0`
path=`realpath $path`
JELLY=$path"/jellyfish-linux"
# calculate kmercount-count table and ordered print it for maternal mers
$JELLY histo -o maternal.histo maternal_mer_counts.jf
# calculate kmercount-count table and ordered print it for paternal mers
$JELLY histo -o paternal.histo pternal_mer_counts.jf

awk -f $path"/find_bounds.awk" maternal.histo > maternal.bounds.txt
awk -f $path"/find_bounds.awk" paternal.histo > paternal.bounds.txt

