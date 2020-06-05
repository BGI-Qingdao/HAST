#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/log/log.h"
#include "common/error/Error.h"

#include "biocommon/fasta/fasta.h"
#include "appcommon/Idx.h"
#include "appcommon/SegmentFa.h"
#include "appcommon/FileName.h"

#include <functional>
#include <map>
#include <string>
#include <sstream>
#include <cassert>

typedef BGIQD::FASTA::Id_Desc_Head FAHEAD;
typedef BGIQD::FASTA::Fasta<FAHEAD> FA ;
typedef BGIQD::FASTA::FastaReader<FA> FA_READER;
typedef BGIQD::APP::Scaff_Seg_Head SEGHEAD ;

int main(int argc , char ** argv)
{
    // -------------------------------------------------------------
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string, fa_1 , "xxx.1.fa");
    DEFINE_ARG_REQUIRED(std::string, fa_2 , "xxx.2.fa");
    DEFINE_ARG_REQUIRED(std::string, idx_1 , "xxx.1.idx");
    DEFINE_ARG_REQUIRED(std::string, idx_2 , "xxx.2.idx");
    DEFINE_ARG_REQUIRED(std::string, prefix , "output prefix\n\
                                                  will output :\n\
                                                      xxx.phb.1.fa\n\
                                                      xxx.phb.2.fa\n\
                                                      xxx.homo.fa");
    END_PARSE_ARGS

    // -------------------------------------------------------------
    BGIQD::APP::FileNames fNames;
    fNames.Init(prefix.to_string());

    // -------------------------------------------------------------
    // functions
    auto load_idx = [](const std::string & line , std::map<unsigned int ,  BGIQD::APP::Idx> & cache) { 
        BGIQD::APP::Idx tmp ;
        tmp.InitFromString(line);
        assert(tmp.is_valid());
        cache[tmp.scaffold_id] = tmp ;
    };
    auto print_seg = [](std::ostream & ost ,
                        const BGIQD::APP::Scaff_Seg_Head & head ,
                         const BGIQD::SEQ::seq & seq ,
                         int start,
                         int end )
    {
        ost<<head.Head()<<'\n';
        BGIQD::SEQ::seq tmp ;
        tmp.AddPartSeq(seq.atcgs.substr(start,end-start));
        ost<<tmp.Seq(60);
    };

    // -------------------------------------------------------------
    // Load xxx.1.idx
    auto in_idx1 =  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(idx_1.to_string());
    if( 0 == in_idx1 ) FATAL("failed to open idx_1 to read");
    std::map<unsigned int , BGIQD::APP::Idx> idx_cache1 ; 
    BGIQD::FILES::FileReaderFactory::EachLine(*in_idx1 , 
            std::bind( load_idx , std::placeholders::_1 , std::ref(idx_cache1))) ;
    delete in_idx1;

    // -------------------------------------------------------------
    // Load xxx.1.fa
    auto in_fa1 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fa_1.to_string());
    if( 0 == in_fa1 ) FATAL("failed to open fa_1 to read");
    FA_READER reader1 ;
    std::vector<FA> buffer1 ;
    reader1.LoadAllFasta(*in_fa1,buffer1);
    delete in_fa1;

    // -------------------------------------------------------------
    // print prefix.phd1.fa
    auto out_1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.phb1_fa());
    if( 0 == out_1 ) FATAL("failed to open xxx.phb.1.fa to write");
    for( const auto & scaff: buffer1 ) {
        unsigned int scaff_id = std::stoul(scaff.head.Id);
        const auto & idx = idx_cache1.at(scaff_id);
        if( ! idx.is_multi() ) continue ;
        auto coords = idx.phase_parts();
        int i = 1 ;
        for( const auto & pair : coords) {
            SEGHEAD tmp ;
            tmp.scaff_id = scaff_id ;
            tmp.seq_index = i ; i+=2 ;
            tmp.phase_id = 1 ;
            print_seg(*out_1,tmp,scaff.seq,pair.first ,pair.second);
        }
    }
    delete out_1;

    // -------------------------------------------------------------
    // print prefix.homo.fa
    auto out_2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.homo_fa());
    if( 0 == out_2 ) FATAL("failed to open xxx.homo.fa to write");
    for( const auto & scaff: buffer1 ) {
        unsigned int scaff_id = std::stoul(scaff.head.Id);
        const auto & idx = idx_cache1.at(scaff_id);
        assert( idx.is_valid()) ;
        auto coords = idx.homo_parts();
        int i = 0 ;
        for( const auto & pair : coords) {
            SEGHEAD tmp ;
            tmp.scaff_id = scaff_id ;
            tmp.seq_index = i ; 
            i+=2 ;
            tmp.phase_id = 0 ;
            print_seg(*out_2,tmp,scaff.seq,pair.first ,pair.second);
        }
    }
    delete out_2;

    // -------------------------------------------------------------
    // release buffers 
    buffer1.clear();
    idx_cache1.clear();

    // -------------------------------------------------------------
    // Load xxx.2.idx
    auto in_idx2 =  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(idx_2.to_string());
    if( 0 == in_idx2 ) FATAL("failed to open idx_2 to read");
    std::map<unsigned int , BGIQD::APP::Idx> idx_cache2; 
    BGIQD::FILES::FileReaderFactory::EachLine(*in_idx2 
            , std::bind( load_idx , std::placeholders::_1 , std::ref(idx_cache2))) ;
    delete in_idx2;

    // -------------------------------------------------------------
    // Load xxx.2.fa 
    auto in_fa2 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fa_2.to_string());
    if( 0 == in_fa2 ) FATAL("failed to open fa_2 to read");
    FA_READER reader2 ;
    std::vector<FA> buffer2 ;
    reader2.LoadAllFasta(*in_fa2,buffer2);
    delete in_fa2;

    // -------------------------------------------------------------
    // print prefix.phd2.fa 
    auto out_3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.phb2_fa());
    if( 0 == out_3 ) FATAL("failed to open xxx.phb.2.fa to write");
    for( const auto & scaff: buffer2 ) {
        unsigned int scaff_id = std::stoul(scaff.head.Id);
        const auto & idx = idx_cache2.at(scaff_id);
        if( ! idx.is_multi() ) continue ;
        auto coords = idx.phase_parts();
        int i = 1 ;
        for( const auto & pair : coords) {
            SEGHEAD tmp ;
            tmp.scaff_id = scaff_id ;
            tmp.seq_index = i ; i+=2 ;
            tmp.phase_id = 2 ;
            print_seg(*out_3,tmp,scaff.seq,pair.first ,pair.second);
        }
    }
    delete out_3;

    // -------------------------------------------------------------
    // exit
    buffer2.clear();
    idx_cache2.clear();
    return 0;
}
