BEGIN{
    used=0;
    total=0;
}
{
    if( FILENAME == ARGV[1] ) {
        s[$1]=1
    } else { 
        if( FNR%4==1 && NF>1 ){
            total = total + 1 ;
            if ( $2 in s ){ 
                print $0 ;
                used = used + 1;
                c=1;
            } else {
                c=0
            }
        } else {
            if(c==1) { 
                print $0 ; 
            }
        }
    }
}
END{
    printf("use %d from %d\n",used,total); >>"filter_reads.log"
}
