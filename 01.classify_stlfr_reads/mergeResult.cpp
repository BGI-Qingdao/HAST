#include <getopt.h>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
float g_hap0_fac =1.0f;
float g_hap1_fac =1.0f;

//
// barcode haplotype relate functions
//
struct BarcodeCache {
    std::map<std::string, std::map<int,int>> barcode_haps;
    void IncrBarcodeHaps(const std::string & barcode , int hap,int incr=1){
        if( barcode_haps[barcode].find(hap) ==  barcode_haps[barcode].end() )
            barcode_haps[barcode][hap] = 0 ;
        barcode_haps[barcode][hap] +=incr ;
    }
    void AddLine(const std::string & line){
        std::string barcode ;
        int type ;
        int hap0 ;
        int hap1 ;
        std::istringstream is(line);
        is>>barcode>>type>>hap0>>hap1;
        IncrBarcodeHaps(barcode,0,hap0);
        IncrBarcodeHaps(barcode,0,hap1);
    }
};

int getHap(const std::string & barcode , const std::map<int,int> & data){
    if( barcode == "0_0_0" || barcode == "0_0" || barcode == "0" )
        return -1;
    if( data.find(0) != data.end() &&  data.find(1) != data.end() ){
        double df0 =  double(data.at(0)) ;/// double(g_kmers[0].size());
        double df1 =  double(data.at(1)) ;/// double(g_kmers[1].size());
        df0 *= g_hap0_fac;
        df1 *= g_hap1_fac;
        if( df0 > df1 ) return 0 ; 
        if( df1 > df0 ) return 1 ;
        return -1 ;
    } else if ( data.find(0) != data.end() ) {
        if ( data.at(0) > 0 ) return 0 ;
        return -1 ;
    } else if ( data.find(1) != data.end() ) {
        if ( data.at(1) > 0 ) return 1 ;
        return -1 ;
    } else {
        return -1 ;
    }
}

int getHapCount(const std::map<int,int> & data , int hap){
    if ( data.find(hap) != data.end() ) return data.at(hap);
    return 0 ;
}

void printBarcodeInfos(const BarcodeCache& g_barcode_haps){
    for(const auto & pair : g_barcode_haps.barcode_haps){
        std::cout<<pair.first;
        const auto & data=pair.second;
        std::cout<<'\t'<<getHap(pair.first,data);
        std::cout<<'\t'<<getHapCount(data,0);
        std::cout<<'\t'<<getHapCount(data,1)<<'\n';
        //std::cout<<'\t'<<getHapCount(data,-1)<<'\n';
    }
}
void printUsage(){
    std::cerr<<"Usage :"<<'\n';
    std::cerr<<"    mergeResult --input i1.txt [--input i2.txt --input ...] [OPTIONS]"<<'\n';
    std::cerr<<"\n";
    std::cerr<<"Options:\n";
    std::cerr<<"    --input         result file to be merged.\n";
    std::cerr<<"    --weight0       weight of hap0\n";
    std::cerr<<"                    notice: weight0=(#hap0's unshared kmers)*(weight0 of classify)\n";
    std::cerr<<"    --weight1       weight of hap1\n";
    std::cerr<<"                    notice: weight1=(#hap1's unshared kmers)*(weight1 of classify)\n";
    return ;
}

int main( int argc , char ** argv){
    static struct option long_options[] = {
        {"input",  required_argument, NULL, 'i'},
        {"weight0" ,required_argument,NULL, 'w'},
        {"weight1" ,required_argument,NULL, 'u'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "p:w:u:h";
    std::vector<std::string> inputs;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'i':
                inputs.push_back(std::string(optarg));
                break;
            case 'u':
                g_hap1_fac=atof(optarg);
                break;
            case 'w':
                g_hap0_fac=atof(optarg);
                break;
            case 'h':
            default :
                printUsage();
                return -1;
        }
    }
    if( inputs.empty()) {
        printUsage();
        return -1;
    }
    BarcodeCache cache;
    for( int i = 0 ; i<(int)inputs.size() ; i++ ){
        std::ifstream in(inputs[i]);
        if( ! in.is_open() ){
            std::cerr<<"ERROR : failed to open "<<inputs[i]<<" to read !!! exit ... "<<std::endl;
            return -1;
        }
        std::string line ;
        while(!std::getline(in,line).eof()){
            cache.AddLine(line);
        }
        in.close();
    }
    printBarcodeInfos(cache);
    return 0;
}

