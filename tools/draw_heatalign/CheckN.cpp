#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

struct NPos{
    std::string seq_name;
    int start ;
    int end ;

    void InitFromString(const std::string & line){
        int det=0;
        for(char c:line) if(c=='\t')det++;
        if(det<2) {
            std::cerr<<"nzone is invalid:"<<std::endl;
            std::cerr<<line<<std::endl;
            std::cerr<<"please use \\t to seperate columns!!!"<<std::endl;
            std::cerr<<"exit ..."<<std::endl;
        }
        std::istringstream ist(line);
        ist>>seq_name>>start>>end;
    }
};

std::map<std::string ,std::vector<NPos> >N_cahces;

struct AlignBlock {
    std::string m_line ;
    int ref_start;
    int ref_end ;
    std::string ref_name;
    std::string query_name;
    int query_start;
    int query_end;
    float idy;
    bool orient;
    int maped_len() const {
        return ref_end - ref_start + 1 ;
    }
    void InitFromString(const std::string & line){
        m_line=line;
        int det=0;
        for(char c:line) if(c=='\t')det++;
        if(det<6) {
            std::cerr<<"align info is invalid:"<<std::endl;
            std::cerr<<line<<std::endl;
            std::cerr<<"please use \\t to seperate columns!!!"<<std::endl;
            std::cerr<<"exit ..."<<std::endl;
        }
        std::istringstream ist(line);
        if( det == 6 ){
            ist>>ref_name>>ref_start>>ref_end>>query_name>>query_start>>query_end>>idy;
            if( query_start < query_end ) orient = true;
            else orient = false ;
        }
        if( det >= 7 ) {
            char o;
            ist>>ref_name>>ref_start>>ref_end>>query_name>>query_start>>query_end>>o>>idy;
            if(o == '+' ) {
                orient = true ;
            } 
            else {
                orient = false ;
                if( query_start < query_end )
                    std::swap(query_start,query_end);
            }
        }
    }
};

std::vector<AlignBlock> block;

int main(int argc , char ** argv){
    if( argc != 3 ) {
        std::cerr<<"Usage : CheckN xxx.align.txt xxx.nzone.txt"<<std::endl;
        return 1;
    }
    std::string align_txt=std::string(argv[1]);
    std::string nzone_txt=std::string(argv[2]);
    return 1 ;
}
