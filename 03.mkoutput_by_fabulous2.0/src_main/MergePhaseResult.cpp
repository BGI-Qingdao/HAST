#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/log/log.h"
#include "common/error/Error.h"

#include "appcommon/FileName.h"
#include "appcommon/SegmentFa.h"

#include <functional>
#include <map>
#include <string>
#include <sstream>
#include <cassert>
#include <set>

std::set<std::string> final_homo_ids;

struct MergeIdsElement {
    enum TrioBinResult {
        TrioBinFather = 1,
        TrioBinMother =2 ,
        TrioBinHomo ,
    } trio_result , trio_paired_result;

    static TrioBinResult Oppo( TrioBinResult t) {
        assert( t != TrioBinResult::TrioBinHomo ) ;
        return ( t == TrioBinResult::TrioBinFather ) ?
            TrioBinResult::TrioBinMother 
            : TrioBinResult::TrioBinFather ;
    }
    enum SupernovaResult {
        Type_1 ,
        Type_2 ,
    } super_result ;

    bool isFinalHomo() const { 
        return trio_paired_result == TrioBinResult::TrioBinHomo ;
    }

    void SetHomeBySupernova( TrioBinResult type_1_eq ) {
        if ( ! isFinalHomo() ) return ;
        if ( super_result == Type_1 )
            trio_paired_result = type_1_eq ;
        else
            trio_paired_result = Oppo(type_1_eq);
    }

    std::string line ;
    float weight;
};

struct MergeIdsPair {
    MergeIdsElement super_1 ;
    MergeIdsElement super_2 ;

    void GenTrioBinPairedResult() {
        if( super_1.trio_result != super_2.trio_result
                && super_1.trio_result != MergeIdsElement::TrioBinResult::TrioBinHomo
                && super_2.trio_result != MergeIdsElement::TrioBinResult::TrioBinHomo){
            super_1.trio_paired_result = super_1.trio_result ;
            super_2.trio_paired_result = super_2.trio_result ;
        } else if ( super_1.trio_result == super_2.trio_result ) {
            super_1.trio_paired_result = MergeIdsElement::TrioBinResult::TrioBinHomo;
            super_2.trio_paired_result = MergeIdsElement::TrioBinResult::TrioBinHomo;
        } else {
            if( super_1.trio_result == MergeIdsElement::TrioBinResult::TrioBinHomo ) {
                super_2.trio_paired_result = super_2.trio_result ;
                super_1.trio_paired_result = MergeIdsElement::Oppo(super_2.trio_result);
            } else {
                assert( super_2.trio_result == MergeIdsElement::TrioBinResult::TrioBinHomo ) ; 
                super_1.trio_paired_result = super_1.trio_result ;
                super_2.trio_paired_result = MergeIdsElement::Oppo(super_1.trio_result);
            }
        }
    }

    bool SetHomeBySupernova(  MergeIdsElement::TrioBinResult type_1_eq ) {
        if ( super_2.trio_paired_result == super_1.trio_paired_result ) {
            assert( super_2.trio_paired_result == MergeIdsElement::TrioBinResult::TrioBinHomo);
            assert( super_2.super_result != super_1.super_result );
            if(  super_1.weight > super_2.weight ) {
                assert( super_1.trio_result != MergeIdsElement::TrioBinResult::TrioBinHomo );
                super_1.trio_paired_result = super_1.trio_result ;
                super_2.trio_paired_result = MergeIdsElement::Oppo(super_1.trio_result);
                return false;
            } else if ( super_1.weight < super_2.weight ) {
                assert( super_2.trio_result != MergeIdsElement::TrioBinResult::TrioBinHomo );
                super_2.trio_paired_result = super_2.trio_result ;
                super_1.trio_paired_result = MergeIdsElement::Oppo(super_2.trio_result);
                return false;
            } else {
                super_1.SetHomeBySupernova(type_1_eq);
                super_2.SetHomeBySupernova(type_1_eq);
                return true;
            }
        }
        return false;
    }

    MergeIdsElement::TrioBinResult VoteSupernovaType1() const {
        return super_1.super_result == MergeIdsElement::SupernovaResult::Type_1 ? 
            super_1.trio_paired_result :
            super_2.trio_paired_result ;
    }
    std::string FatherIds() const {
        assert( ! super_1.isFinalHomo() && ! super_2.isFinalHomo() );
        return super_1.trio_paired_result == MergeIdsElement::TrioBinFather ? 
            super_1.line : super_2.line ; 
    }
    std::string MotherIds() const {
        assert( ! super_1.isFinalHomo() && ! super_2.isFinalHomo() );
        return super_1.trio_paired_result == MergeIdsElement::TrioBinFather ? 
            super_2.line : super_1.line ; 
    }
};
//       scaffold_id  ,     block_id , block_info
std::map<unsigned int , std::map<int , MergeIdsPair> > data;

void GenAllTrioBinPairedResult(){
    for( auto & pair : data ) 
        for( auto & pair2 : pair.second ) 
            pair2.second.GenTrioBinPairedResult();
}
void SetAllHomo(MergeIdsElement::TrioBinResult t){
    for( auto & pair : data ) {
        for( auto & pair2 : pair.second ) {
            if( pair2.second.SetHomeBySupernova(t) ) {
                final_homo_ids.insert(pair2.second.super_1.line);
            }
        }
    }
}

MergeIdsElement::TrioBinResult  CountSupernovaType1(){
    std::map<MergeIdsElement::TrioBinResult , int > counts ;
    counts[MergeIdsElement::TrioBinHomo] = 0 ;
    counts[MergeIdsElement::TrioBinFather ] = 0 ;
    counts[MergeIdsElement::TrioBinMother ] = 0 ;
    for( auto & pair : data )
        for( auto & pair2 : pair.second ){
            auto ret = pair2.second.VoteSupernovaType1();
            counts[ret] ++ ;
        }
    int total = 0 ;
    total += counts[MergeIdsElement::TrioBinMother ];
    total += counts[MergeIdsElement::TrioBinFather ];
    total += counts[MergeIdsElement::TrioBinHomo ];
    float father_fac = float(counts[MergeIdsElement::TrioBinFather ] ) / float(total) ;
    float mother_fac = float(counts[MergeIdsElement::TrioBinMother ] ) / float(total) ;
    float homo_fac = float(counts[MergeIdsElement::TrioBinHomo ] ) / float(total) ;
    std::cerr<<" father_fac "<<father_fac<<std::endl;
    std::cerr<<" mother_fac "<<mother_fac<<std::endl;
    std::cerr<<" homo_fac "<<homo_fac<<std::endl;
    return  father_fac >= mother_fac ? MergeIdsElement::TrioBinFather : MergeIdsElement::TrioBinMother ;
}

void loadIds( const std::string & file_name , MergeIdsElement::TrioBinResult t) {
    auto in =  BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file_name);
    if( 0 == in ) FATAL("failed to open ids to read ");
    std::string line ;
    while( ! std::getline( *in , line ).eof() ) {
        BGIQD::APP::Scaff_Seg_Head tmp ;
        std::istringstream ist(line);
        std::string read_name ;
        float weight ;
        ist>>read_name>>weight;
        tmp.Init(read_name);
        MergeIdsElement elem ;
        elem.trio_result = t ;
        elem.line = read_name ;
        elem.weight = weight ;
        if( tmp.phase_id == 1 ){
            elem.super_result = MergeIdsElement::Type_1 ;
            data[tmp.scaff_id][tmp.seq_index].super_1 = elem ; 
        } else if ( tmp.phase_id == 2 ){
            elem.super_result = MergeIdsElement::Type_2 ;
            data[tmp.scaff_id][tmp.seq_index].super_2 = elem ; 
        } else assert(0) ; 
    }
    delete in;
}

BGIQD::APP::FileNames fNames ;
void PrintFinalIds() {
    auto out1 = BGIQD::FILES::FileWriterFactory
        ::GenerateWriterFromFileName(fNames.merge_father_ids());
    if( 0 == out1 ) FATAL("failed to open xxx.merge.father.ids to write");
    for( auto & pair : data ) 
        for( auto & pair2 : pair.second ) 
            (*out1)<<pair2.second.FatherIds()<<'\n';
    delete out1;
    auto out2 = BGIQD::FILES::FileWriterFactory
        ::GenerateWriterFromFileName(fNames.merge_mother_ids());
    if( 0 == out2 ) FATAL("failed to open xxx.merge.mother.ids to write");
    for( auto & pair : data ) 
        for( auto & pair2 : pair.second ) 
            (*out2)<<pair2.second.MotherIds()<<'\n';
    delete out2;
}

void PrintHomoIds(){
    auto out1 = BGIQD::FILES::FileWriterFactory
        ::GenerateWriterFromFileName(fNames.merge_homo_ids());
    if( 0 == out1 ) FATAL("failed to open xxx.merge.homo.ids to write");
    for( auto & name : final_homo_ids) 
        (*out1)<<name<<'\n';
    delete out1;
}

int main(int argc , char ** argv ) {

    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED( std::string , prefix , "output prefix\n\
                will output: \n\
                xxx.merge.father.ids\n\
                xxx.merge.mother.ids");
    DEFINE_ARG_REQUIRED( std::string , father_ids , "phb12.father_idx from triobin output");
    DEFINE_ARG_REQUIRED( std::string , mother_ids , "phb12.mother_idx from triobin output");
    DEFINE_ARG_REQUIRED( std::string , homo_ids , "phb12.homo_idx from triobin output");
    END_PARSE_ARGS;

    fNames.Init(prefix.to_string());

    loadIds( father_ids.to_string() , MergeIdsElement::TrioBinResult::TrioBinFather) ;
    loadIds( mother_ids.to_string() , MergeIdsElement::TrioBinResult::TrioBinMother) ;
    loadIds( homo_ids.to_string()   , MergeIdsElement::TrioBinResult::TrioBinHomo) ;

    GenAllTrioBinPairedResult() ;
    auto type_1_eq = CountSupernovaType1() ;
    SetAllHomo(type_1_eq);

    PrintFinalIds();
    PrintHomoIds();
    return 0 ;
}
