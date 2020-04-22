#include <iostream>
#include <map>
#include <unordered_set>
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
#include "kmer/kmer.h"
int Kmer::overlap = 0 ;
Kmer Kmer::WORDFILTER ;
void logtime(){
    time_t now = time(0);
    char* dt = ctime(&now);
    std::cerr<<dt<<std::endl;
}
double g_hap0_fac = 1.0f ;
double g_hap1_fac = 1.0f ;
//
// load & cache maternal unique kmer & paternal unique kmer
//
std::unordered_set<Kmer> g_kmers[2];

int g_K=0;
void load_kmers(const std::string & file,int index){
    std::ifstream ifs(file);
    std::string line;
    int total_kmer = 0 ;
    if(index==0){
        std::getline(ifs,line);
        g_K = line.size();
        Kmer::InitFilter(g_K);
        g_kmers[index].insert(Kmer::str2Kmer(BaseStr::str2BaseStr(line)));
        total_kmer++;
    }
    while(!std::getline(ifs,line).eof()){
        g_kmers[index].insert(Kmer::str2Kmer(BaseStr::str2BaseStr(line)));
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
        double df0 =  double(data.at(0)) / double(g_kmers[0].size());
        double df1 =  double(data.at(1)) / double(g_kmers[1].size());
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

#define max_buffer 1024
struct Buffer{
    std::array<std::string,max_buffer> heads;
    std::array<std::string,max_buffer> seqs ;
    void Init() { size = 0 ; }
    int size ;
};

struct MultiThread {
    int t_nums ;
    bool end;
    bool busy;
    void Worker(int index){
        int miss = 0;
        int hit = 0 ;
        Buffer buffer;
        //std::pair<std::string,std::string> job;
        while(true){
            locks[index].lock();
            if( caches[index].empty() ){
                miss ++ ;
                busy = false ;
                locks[index].unlock();
                if(end) { 
                    std::cerr<<"thread="<<index<<" miss="<<miss<<" hit="<<hit<<std::endl;
                    return ;
                }
                std::this_thread::sleep_for(std::chrono::microseconds(10));
                continue;
            }
            if( ! caches[index].empty() ){
                hit ++ ;
                if( caches[index].size() > 300 ) busy = true ;
                else if ( caches[index].size() < 50 ) busy = false ;
                std::swap(buffer,caches[index].top());
                caches[index].pop();
                locks[index].unlock();
                for( int i = 0 ; i < buffer.size ; i ++ )
                    process_reads(buffer.heads.at(i),buffer.seqs.at(i),index);
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
            caches.push_back(std::stack<Buffer>());
            //caches.push_back(std::stack< std::pair<std::string,std::string> >());
            threads[i] = new std::thread([this,i](){int index=i ;Worker(index); });
        }
    }
    ~MultiThread(){
        delete [] threads;
        delete [] locks;
        delete [] barcode_caches;
    }
    void process_reads(const std::string & head ,
                         const std::string & read , int index) {
        int vote[2] ; vote[0]=0;vote[1]=0;
        std::string barcode = parseName(head);
        std::vector<Kmer> kmers=Kmer::chopRead2Kmer(BaseStr::str2BaseStr(read));
        for(int i = 0 ; i <(int)kmers.size();i++){
            if( g_kmers[0].find(kmers.at(i)) != g_kmers[0].end() ){
                vote[0] ++ ;
            }
            if( g_kmers[1].find(kmers.at(i)) != g_kmers[1].end() ){
                vote[1] ++ ;
            }
        }
        if( vote[0] > 0 )
            barcode_caches[index].IncrBarcodeHaps(barcode,0,vote[0]);
        if( vote[1] > 0 )
            barcode_caches[index].IncrBarcodeHaps(barcode,1,vote[1]);
        if( vote[0] == 0 && vote[1] == 0 )
            barcode_caches[index].IncrBarcodeHaps(barcode,-1,1);
    }
    //void submit(const std::string & head ,const std::string & seq){
    void submit(const Buffer & buffer ){
        static long index = 0;
        while( busy ) { std::this_thread::sleep_for(std::chrono::seconds(1));}
        index ++ ;
        int id = index % t_nums;
        locks[id].lock();
        caches[id].push(buffer);//std::make_pair(head,seq));
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
    //std::vector<std::stack< std::pair<std::string,std::string> >>  caches;
    std::vector<std::stack<Buffer> > caches;
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
    Buffer buffer ;
    buffer.Init();
    while(!std::getline(*in,head).eof()){
        std::getline(*in,seq);
        //mt.submit(head,seq);
        buffer.heads.at(buffer.size) = head ;
        buffer.seqs.at(buffer.size) = seq;
        buffer.size ++ ;
        if ( buffer.size >= max_buffer ){
            mt.submit(buffer);
            buffer.Init();
        }
        std::getline(*in,tmp);
        std::getline(*in,tmp);
    }
    if( buffer.size > 0 ){
        mt.submit(buffer);
        buffer.Init();
    }
    mt.end=true;
    mt.wait();
    mt.collectBarcodes();
    data.Add(mt.final_data);
}

void printUsage(){
    std::cerr<<"Uasge :\n\tclassify --hap0 hap0 --hap1 hap1 --read read1.fq [--read read2.fq] [--thread t_num (8 default) ] [--weight0 w0 (1.0 default) ] [--weight1 w1 (1.0 default) ] "<<std::endl;
    std::cerr<<"output format: \n\tbarcode haplotype(0/1/-1) kmer_count_hap0 kmer_count_hap1"<<std::endl;
    std::cerr<<"notice : --read accept file in gzip format , but file must end by \".gz\""<<std::endl;
}
void InitAdaptor(){
    std::string r1("CTGTCTCTTATACACATCTTAGGAAGACAAGCACTGACGACATGATCACCAAGGATCGCCATAGTCCATGCTAAAGGACGTCAGGAAGGGCGATCTCAGG");
    std::string r2("TCTGCTGAGTCGAGAACGTCTCTGTGAGCCAAGGAGTTGCTCTGGCGACGGCCACGAAGCTAACAGCCAATCTGCGTAACAGCCAAACCTGAGATCGCCC");
    std::vector<Kmer> kmers = Kmer::chopRead2Kmer(BaseStr::str2BaseStr(r1));
    for(int i = 0 ; i <(int)kmers.size();i++){
        if( g_kmers[0].find(kmers.at(i)) != g_kmers[0].end() ){
            g_kmers[0].erase(kmers.at(i));
            std::cerr<<" INFO : erase a adaptor kmer from hap 0 ; kmer= "<<BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers.at(i)))<<std::endl;
        }
        if( g_kmers[1].find(kmers.at(i)) != g_kmers[1].end() ){
            g_kmers[1].erase(kmers.at(i));
            std::cerr<<" INFO : erase a adaptor kmer from hap 1 ; kmer= "<<BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers.at(i)))<<std::endl;
        }
    }
    std::vector<Kmer> kmers2 = Kmer::chopRead2Kmer(BaseStr::str2BaseStr(r2));
    for(int i = 0 ; i <(int)kmers2.size();i++){
        if( g_kmers[0].find(kmers2.at(i)) != g_kmers[0].end() ){
            g_kmers[0].erase(kmers2.at(i));
            std::cerr<<" INFO : erase a adaptor kmer from hap 0 ; kmer= "<<BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers2.at(i)))<<std::endl;
        }
        if( g_kmers[1].find(kmers2.at(i)) != g_kmers[1].end() ){
            g_kmers[1].erase(kmers2.at(i));
            std::cerr<<" INFO : erase a adaptor kmer from hap 1 ; kmer= "<<BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers2.at(i)))<<std::endl;
        }
    }
}

void TestAll(){
    assert(parseName("VSDSDS#XXX_xxx_s/1")=="XXX_xxx_s");
    Kmer::InitFilter(5);
    auto str1=BaseStr::str2BaseStr("AGCTC");
    int  t1[] = { '\000','\003','\001','\002','\001'};
    for( int i = 0 ; i <5 ; i++ ) assert(str1[i] == t1[i]);
    auto str2 = BaseStr::str2BaseStr("GAGCT");
    int  t2[] = {'\003','\000','\003','\001','\002'};
    for( int i = 0 ; i <5 ; i++ ) assert(str2[i] == t2[i]);
                         // 00 1101 1001
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("AGCTC")).low == 0xD9);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("AGCTC")).high == 0);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("GAGCT")).low == 0xD9);
    assert(Kmer::str2Kmer(BaseStr::str2BaseStr("GAGCT")).high == 0);
    auto kmers = Kmer::chopRead2Kmer(BaseStr::str2BaseStr("GAGCTA"));
    assert(kmers.size() == 2 );
    // GAGCT->AGCTC 00 1101 1001 -> 0xD9
    // AGCTA        00 1101 1000 -> 0xD8
    assert(kmers[0].high == 0 );
    assert(kmers[0].low == 0xD9 );
    assert(kmers[1].high == 0 );
    assert(kmers[1].low == 0xD8 );
    std::string k1 = BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers[0]));
    std::string k2 = BaseStr::BaseStr2Str(Kmer::ToBaseStr(kmers[1]));
    assert(k1 == "AGCTC");
    assert(k2 == "AGCTA");
}

//
// Main function
//

int main(int argc ,char ** argv ){
    TestAll();
    static struct option long_options[] = {
        {"hap0",  required_argument, NULL, 'p'},
        {"hap1",  required_argument, NULL, 'm'},
        {"read", required_argument,  NULL, 'r'},
        {"thread",required_argument, NULL, 't'},
        {"weight0",required_argument, NULL, 'w'},
        {"weight1",required_argument, NULL, 'u'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "p:m:l:r:t:w:u:h";
    std::string hap0 , hap1 ;
    std::vector<std::string> read;
    int t_num=8;
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
    if( hap0 == "" || hap1 == "" || read.empty()|| t_num< 1) {
        printUsage();
        return -1;
    }
    std::cerr<<"__START__"<<std::endl;
    std::cerr<<" use hap0 weight "<<g_hap0_fac<<std::endl;
    std::cerr<<" use hap1 weight "<<g_hap1_fac<<std::endl;
    logtime();
    std::cerr<<"__load hap0 kmers__"<<std::endl;
    load_kmers(hap0,0);
    std::cerr<<"__load hap1 kmers__"<<std::endl;
    load_kmers(hap1,1);
    InitAdaptor();
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
