#!/bin/bash

path=`dirname $0`

# calculate kmercount-count table and ordered print it for maternal mers
awk -F '>' '{if(NF>1) a[$2]=a[$2]+1;} END{for(item in a){printf("%d %d\n",item,a[item]);}}' maternal.mer.fa| sort -n >maternal.kmercount.count.txt
# calculate kmercount-count table and ordered print it for paternal mers
awk -F '>' '{if(NF>1) a[$2]=a[$2]+1;} END{for(item in a){printf("%d %d\n",item,a[item]);}}' paternal.mer.fa| sort -n >paternal.kmercount.count.txt

awk -f $path"/find_bounds.awk" maternal.kmercount.count.txt > maternal.bounds.txt
awk -f $path"/find_bounds.awk" paternal.kmercount.count.txt > paternal.bounds.txt

