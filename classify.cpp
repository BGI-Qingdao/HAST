#include <iostream>
#include <map>
#include <set>
#include <fstream>
#include <cassert>
#include <ctime>
#include <thread>
#include <stack>
#include <mutex>
#include <vector>
#include <chrono>
#include <getopt.h>
#include "gzstream/gzstream.h"

void logtime(){
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cerr<<dt<<std::endl;
}
//
// load & cache maternal unique kmer & paternal unique kmer
//
std::set<std::string> g_kmers[2];
int g_K=0;
void load_kmers(const std::string & file,int index){
    std::ifstream ifs(file);
    std::string line;
    int total_kmer = 0 ;
    if(index==0){
        std::getline(ifs,line);
        g_K = line.size();
        g_kmers[index].insert(line);
        total_kmer++;
    }
    while(!std::getline(ifs,line).eof()){
        g_kmers[index].insert(line);
        total_kmer++;
    }
    std::cerr<<"Recorded "<<total_kmer<<" haplotype "<<index<<" specific "<<g_K<<"-mers\n"; 
}
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
    void Add(const BarcodeCache & other){
        for(const auto & pair : other.barcode_haps ) {
            for( const auto & pair1 : pair.second ){
                IncrBarcodeHaps(pair.first,pair1.first,pair1.second);
            }
        }
    }
};

int getHap(const std::string & barcode , const std::map<int,int> & data){
    if( barcode == "0_0_0" || barcode == "0_0" || barcode == "0" )
        return -1;
    if( data.find(0) != data.end() &&  data.find(1) != data.end() ){
        if( data.at(0) > data.at(1) ) return 0 ; 
        if( data.at(1) > data.at(0) ) return 1 ;
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
        std::cout<<'\t'<<getHapCount(data,1);
        std::cout<<'\t'<<getHapCount(data,-1)<<'\n';
    }
}

//
// kmer relate functions
//
std::map<char,char> g_oppo ;
void InitMap(){
    g_oppo['a']='T';
    g_oppo['A']='T';
    g_oppo['g']='C';
    g_oppo['G']='C';
    g_oppo['c']='G';
    g_oppo['C']='G';
    g_oppo['t']='A';
    g_oppo['T']='A';
}
std::string reverse_complement(const std::string & kmer){
    std::string ret = kmer;
    for(int i = 0 ; i< (int)kmer.size(); i++)
        ret[kmer.size()-i-1] = g_oppo[kmer[i]];
    return ret ;
}

std::string get_cannonical(const std::string & kmer){
   std::string rc = reverse_complement(kmer);
   if (rc < kmer)
      return rc;
   else
      return kmer;
}

//
//reads relate functions
//

// @return : barcode
// A stLFR read's head looks like :
//   @V300017823L1C001R051096800#203_1533_1069/1
//                barcode str :  203_1533_1069
std::string parseName(const std::string & head){
    int s=-1, e=-1;
    for( int i = 0 ; i< (int)head.size(); i++ ){
        if( head[i] == '#' ) s=i;
        if( head[i] == '/' ) e=i;
    }
    return head.substr(s+1,e-s-1);
}

struct MultiThread {
    int t_nums ;
    bool end;
    bool busy;
    void Worker(int index){
        std::pair<std::string,std::string> job;
        while(true){
            locks[index].lock();
            if( caches[index].empty() ){
                busy = false ;
                locks[index].unlock();
                if(end) return ;
                std::this_thread::sleep_for(std::chrono::microseconds(10));
                continue;
            }
            if( ! caches[index].empty() ){
                job=caches[index].top();
                if( caches[index].size() > 10000000 ) busy = true ;
                else if ( caches[index].size() < 300000 ) busy = false ;
                caches[index].pop();
                locks[index].unlock();
                process_reads(job.first,job.second,index);
            } else 
                locks[index].unlock();
        }
    }
    MultiThread(int t_num){
        t_nums = t_num ;
        barcode_caches = new BarcodeCache[t_num];
        locks = new std::mutex[t_num];
        threads = new std::thread*[t_num];
        busy = false;
        end=false;
        for(int i = 0 ; i< t_num ; i++){
            caches.push_back(std::stack< std::pair<std::string,std::string> >());
            threads[i] = new std::thread([this,i](){int index=i ;Worker(index); });
        }
    }
    ~MultiThread(){
        delete [] threads;
        delete [] locks;
        delete [] barcode_caches;
    }
    void process_reads(const std::string & head ,
                         const std::string & seq , int index) {
        int vote[2] ; vote[0]=0;vote[1]=0;
        for(int i = 0 ; i <(int)seq.size()-g_K+1;i++){
            std::string kmer = get_cannonical(seq.substr(i,g_K));
            if( g_kmers[0].find(kmer) != g_kmers[0].end() )
                vote[0] ++ ;
            if( g_kmers[1].find(kmer) != g_kmers[1].end() )
                vote[1] ++ ;
        }
        std::string barcode = parseName(head);
        if( vote[0] > 0 && vote[0] > vote[1] )
            barcode_caches[index].IncrBarcodeHaps(barcode,0);
        else if( vote[1] > 0 && vote[1] > vote[0] )
            barcode_caches[index].IncrBarcodeHaps(barcode,1);
        else
            barcode_caches[index].IncrBarcodeHaps(barcode,-1);
    }
    void submit(const std::string & head ,const std::string & seq){
        static long index = 0;
        while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
        index ++ ;
        int id = index % t_nums;
        locks[id].lock();
        caches[id].push(std::make_pair(head,seq));
        locks[id].unlock();
    }
    void wait(){
        for(int i = 0 ; i <t_nums; i++){
            threads[i]->join();
            delete threads[i];
        }
    }
    void collectBarcodes(){
        for(int i = 0 ; i<t_nums ;i++)
            final_data.Add(barcode_caches[i]);
    }
    std::vector<std::stack< std::pair<std::string,std::string> >>  caches;
    std::mutex * locks;
    std::thread ** threads; 
    BarcodeCache * barcode_caches;
    BarcodeCache final_data;
};

void processFastq(const std::string & file,int t_num,BarcodeCache& data){
    std::string head;
    std::string seq;
    std::string tmp;
    MultiThread mt(t_num);
    std::istream *in ;
    bool gz_file = false;
    if( file.size() > 3 ) {
        int end=file.size() ;
        if( file[end-3] == '.' && file[end-2] == 'g' && file[end-1]=='z' ) {
            gz_file = true ;
        }
    }
    if ( gz_file )
        in = new igzstream(file.c_str());
    else 
        in = new std::ifstream(file);
    while(!std::getline(*in,head).eof()){
        std::getline(*in,seq);
        mt.submit(head,seq);
        std::getline(*in,tmp);
        std::getline(*in,tmp);
    }
    mt.end=true;
    mt.wait();
    mt.collectBarcodes();
    data.Add(mt.final_data);
}

void printUsage(){
    std::cerr<<"Uasge :\n\tclassify --hap0 hap0 --hap1 hap1 --read read1.fq [--read read2.fq] [--thread t_num]"<<std::endl;
    std::cerr<<"output format: \n\tbarcode haplotype(0/1/-1) read_count_hap0 read_count_hap1 read_count_hap-1"<<std::endl;
    std::cerr<<"notice : --read accept file in gzip format , but file must end by \".gz\""<<std::endl;
}
//
// Main function
//

int main(int argc ,char ** argv ){
    static struct option long_options[] = {
        {"hap0",  required_argument, NULL, 'p'},
        {"hap1",  required_argument, NULL, 'm'},
        {"read", required_argument,  NULL, 'r'},
        {"thread",required_argument, NULL, 't'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "p:m:l:r:t:h";
    std::string hap0 , hap1 ;
    std::vector<std::string> read;
    int t_num=1;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'p':
                hap0 = std::string(optarg);
                break;
            case 'm':
                hap1 = std::string(optarg);
                break;
            case 'r':
                read.push_back(std::string(optarg));
                break;
            case 't':
                t_num = atoi(optarg);
                break;
            case 'h':
            default :
                printUsage();
                return -1;
        }
    }
    if( hap0 == "" || hap1 == "" || read.empty()|| t_num< 1) {
        printUsage();
        return -1;
    }
    InitMap();
    assert(parseName("VSDSDS#XXX_xxx_s/1")=="XXX_xxx_s");
    assert(get_cannonical("AGCTA")=="AGCTA");
    assert(get_cannonical("TGCTT")=="AAGCA");
    std::cerr<<"__START__"<<std::endl;
    logtime();
    std::cerr<<"__load hap0 kmers__"<<std::endl;
    load_kmers(hap0,0);
    std::cerr<<"__load hap1 kmers__"<<std::endl;
    load_kmers(hap1,1);
    logtime();
    BarcodeCache data;
    for(const auto r : read ){
        std::cerr<<"__process read: "<<r<<std::endl;
        processFastq(r,t_num,data);
        logtime();
        std::cerr<<"__process read done__"<<std::endl;
    }
    std::cerr<<"__print result__"<<std::endl;
    printBarcodeInfos(data);
    logtime();
    std::cerr<<"__END__"<<std::endl;
}
