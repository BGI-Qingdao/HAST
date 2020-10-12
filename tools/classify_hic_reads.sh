#!/bin/bash


if [[ $# != 2 ]] ; then
    echo "usage: ./classify_hic_reads.sh <r1.sam> <r2.sam>"
    exit
fi

R1SAM=$1
R2SAM=$2

function get_infos(){
    printf "read-name\tflag\tidy\texact-match-len\ttotal-match-len\n"

    #'if(/NM:i:(\d+)/){
    #     $n=$1;
    #     $m=$g=$o=0;
    #     $m+=$1 while/(\d+)M/g;
    #     $g+=$1,++$o while/(\d+)[ID]/g;
    #     chomp ;
    #     @a=split;
    #     print($a[0],"\t",$a[1],"\t",1-($n-$g+$o)/($m+$o) ,"\t",($m+$o)-($n-$g+$o),  "\t",($m+$o),"\n")
    #}
    #else{
    #     print($a[0],"\t0\t0\t0\t0\n")
    #}'
    grep -v "^@" $1 | cut -f 1,2,6,12 | perl -ane  'if(/NM:i:(\d+)/){ $n=$1;$m=$g=$o=0;$m+=$1 while/(\d+)M/g;$g+=$1,++$o while/(\d+)[ID]/g; chomp ;@a=split;print($a[0],"\t",$a[1],"\t",1-($n-$g+$o)/($m+$o) ,"\t",($m+$o)-($n-$g+$o),  "\t",($m+$o),"\n")} else{print($a[0],"\t0\t0\t0\t0\n")} '
}

function get_scores(){
    #
    # awk
    #'BEGIN{
    #       name="";score=0;
    #}
    #{
    #    if(NR>1){
    #        n=$1;
    #        if(n!=name && name!="")
    #        {
    #            printf("%s\t%f\n",name,score);
    #            score=0;
    #        }
    #        name=n;
    #        if($2<256 && $2>0)
    #        {
    #           score+=( 3*(log($3)/log(10)) + (log($5)/log(10)) );
    #        }
    #    }
    #}' $1

    awk 'BEGIN{name="";score=0;}{ if(NR>1){ n=$1; if(n!=name && name!=""){printf("%s\t%f\n",name,score);score=0;} name=n;if($2<256&& $2>0){score+=(3*(log($3)/log(10))+(log($5)/log(10)) );} } }' $1
}

get_infos $R1SAM >r1.info
get_infos $R2SAM >r2.info

get_scores r1.info >r1.scores
get_scores r2.info >r2.scores

#sort -k 1b,1 r1.scores >r1.sort.scores
#sort -k 1b,1 r2.scores >r2.sort.scores
#join -o 1.1 2.1 1.2  2.2  -e 0 -a1 -a2 r1.sort.scores r2.sort.scores  >read.scores 2>join.err
join -o 1.1 2.1 1.2  2.2  -e 0 -a1 -a2 r1.scores r2.scores  >read.scores 2>join.err
#{
#   if($1==0)
#        name =$2;
#   else
#       name = $1 ;
#   if($3>$4)
#       print name>"paternal.reads";
#   else if ( $4>$3)
#       print name>"maternal.reads" ;
#   else
#       print name>"homo.reads";
#}
awk '{if($1==0) name =$2; else name = $1 ; if($3>$4)print name>"paternal.reads"; else if ( $4>$3) print name>"maternal.reads" ; else print name>"homo.reads"}' read.scores

