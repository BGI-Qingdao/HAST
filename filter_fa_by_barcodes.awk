{
    if( FILENAME == ARGV[1] ) {
        s[$1]=1
    } else { 
        if( FNR%2==1 && NF>1 ){
            if ( $2 in s ){ 
                print $0 ;
                c=1;
            } else {
                c=0
            }
        } else {
            if(c==1) { 
                print $0 ; 
                c=0;
            }
        }
    }
}
