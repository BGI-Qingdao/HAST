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

///////////////////////////////////////////////////////////
// typedef & structures
///////////////////////////////////////////////////////////
typedef BGIQD::APP::Scaff_Seg_Head SEGHEAD ;
typedef BGIQD::FASTA::Fasta<SEGHEAD> FA ;
typedef BGIQD::FASTA::FastaReader<FA> FA_READER;

enum PhasedResult {
    Father,
    Mother,
} phased_result ;

struct IBlock {
    virtual const BGIQD::SEQ::seq *getSeq(PhasedResult type) const = 0 ;
    virtual ~IBlock() {} ;
};

struct HomoBlock : public IBlock {
    FA * block ;
    virtual const BGIQD::SEQ::seq *getSeq(PhasedResult /**type**/ ) const final {
        return &(block->seq);
    }
    virtual ~HomoBlock() final {}
};

struct PhasedBlock : public IBlock {
    FA * father_block;
    FA * mother_block;
    virtual const BGIQD::SEQ::seq * getSeq(PhasedResult type ) const final {
        return type == PhasedResult::Father ? &(father_block->seq) : &(mother_block->seq);
    }
    virtual ~PhasedBlock() final {}
};

struct NewScaff {
    unsigned int id ;
    HomoBlock * getHomePtr(int block_id){
        assert( block_id %2 == 0 );
        if( blocks.find(block_id) == blocks.end() ) {
            HomoBlock* tmp = new HomoBlock() ;
            blocks[block_id] = tmp ;
            return tmp ;
        } else {
            return dynamic_cast<HomoBlock*>(blocks.at(block_id)); 
        }
    }

    PhasedBlock * getPhasedPtr(int block_id ){
        assert( block_id %2 == 1 );
        if( blocks.find(block_id) == blocks.end() ) {
            PhasedBlock * tmp = new PhasedBlock() ;
            blocks[block_id] = tmp ;
            return tmp ;
        } else {
            return dynamic_cast<PhasedBlock*>(blocks.at(block_id)); 
        }
    }
    void getSeq( PhasedResult type , std::vector<int> & idx, BGIQD::SEQ::seq & ret ) const{
        assert(ret.atcgs.empty());
        assert(idx.empty());
        idx.push_back(0);
        assert(blocks.size() %2 == 1 );
        for( int i = 0 ; i < (int)blocks.size() ; i ++ ) {
            try {
                const IBlock * block_ptr = blocks.at(i) ;
                ret.AddPartSeq(block_ptr->getSeq(type)->atcgs);
                idx.push_back(ret.Len());
            }catch(...){
                std::cerr<<"some blocks was missing for scaff_id = "<<id<<" and block_id = "<<i<<std::endl;
                FATAL("block not exsit !!!");
            }
        }
    }
    ~NewScaff() {
        for( const auto & pair : blocks ) 
            delete pair.second ;
    }
    private:
        //   block_id , data_ptr
        std::map< int , IBlock * > blocks;
};

///////////////////////////////////////////////////////////
// global variables
///////////////////////////////////////////////////////////
BGIQD::APP::FileNames fNames;
//       scaff_id,      block_id,      phase_id ,fa
std::map<unsigned int , std::map<int , std::map< int ,FA > > > phb_chace ;
std::map<unsigned int , NewScaff> scaffs ;
std::map<SEGHEAD,BGIQD::SEQ::seq *> supplement_seqs;
std::string prefered_name ;
///////////////////////////////////////////////////////////
// functions 
///////////////////////////////////////////////////////////
void LoadFas(){
    auto load_fa = [&](std::istream & ist) ->void {
        FA_READER reader ;
        std::vector<FA> buffer ;
        reader.LoadAllFasta(ist,buffer);
        for( const auto & fa : buffer)
            phb_chace[fa.head.scaff_id][fa.head.seq_index][fa.head.phase_id] = fa ;
        buffer.clear() ;
    };
    // load phb.1.fa
    auto in1=  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.phb1_fa());
    if( 0 == in1) FATAL("failed to open xxx.phb.1.fa to read");
    load_fa(*in1);
    delete in1 ;
    // load phb.2.fa
    auto in2=  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.phb2_fa());
    if( 0 == in2) FATAL("failed to open xxx.phb.2.fa to read");
    load_fa(*in2);
    delete in2 ;
    // load homo.fa
    auto in3=  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.homo_fa());
    if( 0 == in3) FATAL("failed to open xxx.homo.fa to read");
    load_fa(*in3);
    delete in3 ;
}

void BuildHomeBlocs() {
    for( auto &scaff_pair : phb_chace ) { 
        unsigned int scaff_id = scaff_pair.first ;
        for( auto & block_pair : scaff_pair.second ) {
            int block_id = block_pair.first ;
            for( auto & phase_pair : block_pair.second) {
                int phased_id = phase_pair.first ;
                if( phased_id == 0 ) {
                    scaffs[scaff_id].getHomePtr(block_id)->block = & (phase_pair.second);
                }
            }
        }
    }
}

void BuildPhasedBlocks() {
    auto load_ids = [&](std::istream & ist , PhasedResult type ) -> void {
        std::string line ;
        while( ! std::getline(ist , line ).eof() ) {
            BGIQD::APP::Scaff_Seg_Head tmp ;
            tmp.Init(line);
            try {
                FA * faptr= &(phb_chace.at(tmp.scaff_id).at(tmp.seq_index).at(tmp.phase_id));
                PhasedBlock * pptr = scaffs.at(tmp.scaff_id).getPhasedPtr(tmp.seq_index);
                if( type == PhasedResult::Father )
                    pptr->father_block = faptr;
                else 
                    pptr->mother_block = faptr;
            } catch (...) {
                std::cerr<<"some blocks was missing for scaff_id = "<<tmp.scaff_id<<" and block_id = "<<tmp.seq_index<<std::endl;
                FATAL("block not exsit !!!");
            }
        }
    };
    // load xxx.merge.father.ids
    auto in1 =  BGIQD::FILES::FileReaderFactory
        ::GenerateReaderFromFileName(fNames.merge_father_ids());
    if( 0 == in1 ) FATAL("failed to open xxx.merged.father.ids to read ");
    load_ids(*in1,PhasedResult::Father);
    delete in1;

    // load xxx.merge.mother.ids
    auto in2 =  BGIQD::FILES::FileReaderFactory
        ::GenerateReaderFromFileName(fNames.merge_mother_ids());
    if( 0 == in2 ) FATAL("failed to open xxx.merged.mother.ids to read ");
    load_ids(*in2,PhasedResult::Mother);
    delete in2;
}

void PrintFinalDatas() {
    std::map<unsigned int , std::vector<int>> idx_chace ;
    if( prefered_name == "pat" ) {
        auto out1 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.father_fa());
        if( 0 == out1 ) FATAL("failed to open xxx.father.fa to write");
        for( const auto & pair : scaffs ) {
            (*out1)<<'>'<<pair.first<<'\n';
            std::vector<int> idx;
            BGIQD::SEQ::seq seq ;
            pair.second.getSeq(PhasedResult::Father,idx,seq);
            (*out1)<<seq.Seq(80);
            idx_chace[pair.first]=idx;
        }
        delete out1;
        auto out2 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.father_idx());
        if( 0 == out2 ) FATAL("failed to open xxx.father.idx to write");
        for( const auto & idx : idx_chace ) {
            (*out2)<<idx.first ;
            for( int i : idx.second) (*out2)<<' '<<i;
            (*out2)<<'\n';
        }
        delete out2;
        idx_chace.clear();
    } else {
        auto out3 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.mother_fa());
        if( 0 == out3 ) FATAL("failed to open xxx.mother.fa to write");
        for( const auto & pair : scaffs ) {
            (*out3)<<'>'<<pair.first<<'\n';
            std::vector<int> idx;
            BGIQD::SEQ::seq seq ;
            pair.second.getSeq(PhasedResult::Mother,idx,seq);
            (*out3)<<seq.Seq(80);
            idx_chace[pair.first]=idx;
        }
        delete out3;
        auto out4 = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.mother_idx());
        if( 0 == out4 ) FATAL("failed to open xxx.mother.idx to write");
        for( const auto & idx : idx_chace ) {
            (*out4)<<idx.first ;
            for( int i : idx.second) (*out4)<<' '<<i;
            (*out4)<<'\n';
        }
        delete out4;
    }
}

void PrintFinalUnphasedData(){
    auto out4 = BGIQD::FILES::FileWriterFactory
        ::GenerateWriterFromFileName(fNames.supplement_fa());
    if( 0 == out4 ) FATAL("failed to open xxx.supplement.fa to write");
    for( const auto & pair : supplement_seqs){
        const auto & segname = pair.first;
        const auto & seq = *pair.second ;
        (*out4)<<">scaff_"<<segname.scaff_id<<"_segment_"<<segname.seq_index<<'\n';
        (*out4)<<seq.Seq(80);
    }
    delete out4;
}

void BuildUnphasedBlocks(){
    // load xxx.merge.homo.ids
    auto in1 =  BGIQD::FILES::FileReaderFactory
        ::GenerateReaderFromFileName(fNames.merge_homo_ids());
    if( 0 == in1 ) FATAL("failed to open xxx.merged.homo.ids to read ");
    std::string line ;
    while( ! std::getline((*in1) , line ).eof() ) {
        BGIQD::APP::Scaff_Seg_Head tmp ;
        tmp.Init(line);
        try {
            PhasedBlock * pptr = scaffs.at(tmp.scaff_id).getPhasedPtr(tmp.seq_index);
            if( prefered_name == "pat" )
                supplement_seqs[tmp] = & pptr->mother_block->seq;
            else
                supplement_seqs[tmp] = & pptr->father_block->seq;
        } catch (...) {
            std::cerr<<"some blocks was missing for scaff_id = "<<tmp.scaff_id<<" and block_id = "<<tmp.seq_index<<std::endl;
            FATAL("block not exsit !!!");
        }
    }
    delete in1;
}

///////////////////////////////////////////////////////////
// main
//////////////////////////////////////////////////////////
int main(int argc , char ** argv) {
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED( std::string , prefix , "prefix\n\
                                                        need input : \n\
                                                            xxx.homo.fa\n\
                                                            xxx.phb.1.fa\n\
                                                            xxx.phb.2.fa\n\
                                                            xxx.merge.father.ids\n\
                                                            xxx.merge.mother.ids\n\
                                                            xxx.merge.homo.ids\n\
                                                        will output: \n\
                                                            xxx.father.fa\n\
                                                            xxx.mother.fa\n\
                                                            xxx.father.idx\n\
                                                            xxx.mother.idx");
    DEFINE_ARG_REQUIRED(std::string, prefer , "the prefered branch name. pat or mat");
    END_PARSE_ARGS;

    fNames.Init(prefix.to_string()) ;
    prefered_name = prefer.to_string();
    if( prefered_name != "pat" && prefered_name != "mat" ){
        std::cerr<<"invalid argument! prefer is not pat && not mat! exit ..."<<std::endl;
        exit(1);
    }
    LoadFas();
    BuildHomeBlocs();
    BuildPhasedBlocks();
    PrintFinalDatas();
    BuildUnphasedBlocks();
    PrintFinalUnphasedData();
    return 0 ;
}
