#!/bin/bash

SCRIPT_PATH=`dirname $0`
SCRIPT_PATH=`realpath $SCRIPT_PATH`

perl $SCRIPT_PATH"/split.pl" $*
tags=`seq -w 1 100`
for tag in $tags
do
    file="reads-"$tag".fasta.gz"
    if [[ ! -e $file ]] ; then
        break
    fi
    echo "count $file now ..."
    $SCRIPT_PATH"/meryl" k=21 threads=$CPU memory=$MEMORY \
        count output "./reads-"$tag".meryl" $file || exit 1
done

echo "Merge start ...."
$SCRIPT_PATH"/meryl"  threads=$CPU memory=$MEMORY \
        union-sum output ./reads.final.meryl ./reads-???.meryl

echo "Merge end  ....  ALL DONE "
