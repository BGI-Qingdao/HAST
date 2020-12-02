#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>

struct AGene{
    std::string gene_name;
    int start ;
    int end ;
    std::string m_line ;
    void InitFromString(const std::string & line){
        m_line = line ;
        int det=0;
        for(char c:line) if(c=='\t')det++;
        if(det<2) {
            std::cerr<<"gene.txt is invalid:"<<std::endl;
            std::cerr<<line<<std::endl;
            std::cerr<<"please use \\t to seperate columns!!!"<<std::endl;
            std::cerr<<"exit ..."<<std::endl;
        }
        std::istringstream ist(line);
        ist>>start>>end>>gene_name;
    }
};

std::vector<AGene> GeneCaches;

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

std::vector<AlignBlock> block_caches;

float GeneCov(const AGene & gene ) {
    int total_length = gene.end- gene.start + 1 ;
    int cov_length = 0 ;
    for( const auto & block : block_caches ){
        if( block.ref_start > gene.end || block.ref_end < gene.start ) 
            continue ;
        int cov_s = block.ref_start < gene.start ? gene.start : block.ref_start;
        int cov_e = block.ref_end > gene.end ? gene.end : block.ref_end;
        cov_length += cov_e -cov_s + 1 ;
    }
    return float(cov_length) / float(total_length) ;
}

int main(int argc , char ** argv) {
    if( argc != 3 ) {
        std::cerr<<"Usage : CheckGene xxx.align.txt xxx.genes.txt"<<std::endl;
        return 1;
    }
    std::string align_txt=std::string(argv[1]);
    std::string gene_txt=std::string(argv[2]);
    std::ifstream if1(align_txt);
    std::ifstream if2(gene_txt);

    if ( ! if1.is_open() || ! if2.is_open() ) {
        std::cerr<<"open file failed ! exit ..."<<std::endl;
        return 1 ;
    }
    std::string line ;
    while ( ! std::getline( if2 , line).eof() ) {
        AGene temp;
        temp.InitFromString(line);
        GeneCaches.push_back(temp);
    }
    if2.close();
    while( ! std::getline(if1,line).eof() ) {
        AlignBlock temp;
        temp.InitFromString(line);
        block_caches.push_back(temp);
    }
    if1.close();
    std::set<std::string> gene_names;
    for( const auto & gene : GeneCaches ) {
        if( gene_names.find(gene.gene_name) != gene_names.end() )
            continue ;
        std::cout<<gene.m_line<<'\t'<<GeneCov(gene)<<'\n';
        gene_names.insert(gene.gene_name);
    }
    return 0 ;
}
