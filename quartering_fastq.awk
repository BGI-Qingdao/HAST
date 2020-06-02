BEGIN{
    pa_reads=0;
    ma_reads=0;
    ho_reads=0;
    no_reads=0;
    un_reads=0;
    total=0;
    read_type == -1 ; # 0 : no-barcode ; 1 : paternal  ;
                      # 2 : maternal   ; 3 : homozyote ;
}
{
    if( FILENAME == ARGV[1] ) {
        paternal_barcodes[$1]=1;
    } else if ( FILENAME == ARGV[2] ) {
        maternal_barcodes[$1]=1;
    } else if ( FILENAME == ARGV[3] ) {
        homozygote_barcodes[$1]=1;
    } else {
        if( FNR == 1 ) {
            print FILENAME >>"filter_reads.log";
        }
        if( FNR%4 == 1 ) { # head line of fastq
            total = total + 1 ;
            if ( NF > 1 && $2 != "0_0_0" ) {
                if ( $2 in paternal_barcodes ){
                    pa_reads += 1 ;
                    read_type = 1 ;
                } else if ( $2 in maternal_barcodes ) {
                    ma_reads += 1 ;
                    read_type = 2 ;
                } else if ( $2 in homozygote_barcodes ) {
                    ho_reads += 1 ;
                    read_type = 3 ;
                } else {
                    printf("ERROR : unclassify barcode : %s\n",$2) > "/dev/stderr" ;
                    un_reads+=1 ;
                    read_type = -1 ;
                }
            } else {
                no_reads += 1 ;
                read_type = 0 ;
            }
        }
        if( read_type == 0 ) {
            print $0 >prefix".nobarcode.fastq" ;
        } else if ( read_type == 1 ) {
            print $0 >prefix".paternal.fastq" ;
        } else if ( read_type == 2 ) {
            print $0 >prefix".maternal.fastq" ;
        } else if ( read_type == 3 ) {
            print $0 >prefix".homozygote.fastq" ;
        }
    }
}
END{
    printf("#Total reads                : %d \n",total) >>"filter_reads.log";
    printf("#Reads without barcode      : %d \n",no_reads) >>"filter_reads.log";
    printf("#Paternal reads             : %d \n",pa_reads) >>"filter_reads.log";
    printf("#Maternal reads             : %d \n",ma_reads) >>"filter_reads.log";
    printf("#Homozygote reads           : %d \n",ho_reads) >>"filter_reads.log";
}
