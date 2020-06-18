#!/bin/bash

path=`dirname $0`
path=`realpath $path`
MERYL=$path"/meryl-1.0/meryl"
# calculate kmercount-count table and ordered print it for maternal mers
cd maternal_meryl 
$MERYL histogram  reads.final.meryl >../maternal.histo
cd ..
# calculate kmercount-count table and ordered print it for paternal mers
cd paternal_meryl 
$MERYL histogram  reads.final.meryl >../paternal.histo
cd ..
awk -f $path"/find_bounds.awk" maternal.histo > maternal.bounds.txt
awk -f $path"/find_bounds.awk" paternal.histo > paternal.bounds.txt
