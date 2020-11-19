#include <iostream>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <vector>

///////////////////////////////////////////////////////////
//
// Licence : GPL
// 
// Author  : Lidong Guo
//
// E-mail  : guolidng@genomics.cn
//
///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////
//
// Utils
//
///////////////////////////////////////////////////////////

void PrintUsage(){
std::cerr<<"\n\
 Brief   : Filter all phased SNP from a VCF.\n\
\n\
 Usage   : \n\
   ./PhasedSNP  XXX.vcf  >phased_snp.txt\n\
\n\
 Output  :\n\
   CHROM\tPOS\tN1\tN2\tPS\n\
\n\
    N1 refer to ALT of haplotype1\n\
    N2 refer to ALT of haplotype2\n\
\n\
";
}

/////////////////////////////////////////////////////////

std::vector<std::string>  split( const std::string & str , const char& spliter ) 
{
    std::vector<std::string> ret;
    std::string::size_type pos1 = 0;
    std::string::size_type pos2 = str.find(spliter);
    while( pos2 != str.npos )
    {
        auto item = str.substr( pos1 , pos2-pos1 ) ;
        if( ! item.empty() )
            ret.push_back( item );
        pos1 = pos2 + 1 ;
        pos2 = str.find( spliter , pos1 ) ;
    }
    if( pos1 != str.length() )
    {
        ret.push_back(str.substr(pos1));
    }
    return ret ;
}

///////////////////////////////////////////////////////////
//
// Struct & Global variables
//
///////////////////////////////////////////////////////////
enum V_HType{
    type_0_1 = 0 ,
    type_1_1 = 1 ,
    type_1_2 = 2 ,
};

std::string V_HType_ToStr(V_HType t){
    if(t==V_HType::type_0_1) return "0_1";
    if(t==V_HType::type_1_1) return "1_1";
    if(t==V_HType::type_1_2) return "1_2";
    assert(0);
    return "unknow";
}

enum V_Type {
    isSNP = 0,
    isInDel = 1 ,
    isSV = 2
};

std::string V_Type_ToStr(V_Type t){
    if( t == V_Type::isSNP ) return "SNP" ;
    if( t == V_Type::isInDel ) return "InDel" ;
    if( t == V_Type::isSV ) return "SV" ;
    assert(0);
    return "unknow";
}

// Variant Information
struct VI {
    std::string ref_name;
    int pos ;
    std::string ref;
    std::string alt;
    std::set<std::string> seqs;
    std::string  gt_str;
    std::string  phased_id;
    std::string alt1;
    std::string alt2;
    V_HType htype;
    V_Type type ;
    // Func1
    void InitFromVCF( const std::string & line){
        auto items = split(line,'\t');
        ref_name=items[0];
        pos=std::stoi(items[1]);
        ref=items[3];
        alt=items[4];
        auto v_alts=split(items[4],',');
        for( auto x : v_alts ){
            seqs.insert(x);
        }
        auto datas=split(items[9],':');
        gt_str = datas[0];
        phased_id= datas[1];
        if( datas[0] == "0|1" || datas[0] == "0/1" ||  datas[0] == "1|0" ||  datas[0] == "1/0" ) {
            htype =V_HType::type_0_1;
        }
        else if( datas[0] == "1|1" || datas[0] == "1/1" ) {
            htype =V_HType::type_1_1;
        } else {
            htype =V_HType::type_1_2;
        }
        if( htype == V_HType::type_1_1 ) {
            alt1 = v_alts[0];
            alt2 = v_alts[0];
        } else if ( htype == V_HType::type_0_1 ) {
            if( gt_str == "0/1" || gt_str == "0|1" ){
                alt1 = ref ;
                alt2 = v_alts[0];
            } else {
                alt2 = ref ;
                alt1 = v_alts[0];
            }
        } else {
            if( gt_str == "1/2" || gt_str == "1|2" ){
                alt1 = v_alts[0];
                alt2 = v_alts[1];
            } else {
                alt2 = v_alts[1];
                alt1 = v_alts[0];
            }
        }

        if( htype ==V_HType::type_0_1 ) {
            seqs.insert(ref);
        }
        if(  htype ==V_HType::type_1_1  ) {
            if(seqs.size()!=1)
                std::cerr<<line<<std::endl;
            //notice : HG002 high confidence V4.1 do hit this once! 
            //assert(seqs.size()==1);
        }
        else {
            assert(seqs.size()==2);
        }
        type = V_Type::isSNP;
        if( ref.size() == 1 )  {
            for(const auto & x : seqs ) {
                if(x.size() > 1 ) { type = V_Type::isInDel; }
            }
        } else {
            type = V_Type::isInDel;
        }
        if( type == V_Type::isInDel ) {
            if( ref.size() > 50 )  type =V_Type::isSV ;
            for(const auto & x : seqs ) {
                if(x.size() > 50 )  type =V_Type::isSV ;
            }
        }
    }
};

bool isValidVCF(const std::string & line) {
	return line.size()>0 && line.at(0) != '#' ;
}

int main(int argc , char ** argv)
{
    if(argc != 2 ) { 
        PrintUsage();
        return 0;
    }
    std::string F1_VCF(argv[1]);
    if( F1_VCF == "-h" || F1_VCF == "--help" ){
        PrintUsage();
        return 0;
    }

    std::string line;
    int variant_count = 0;
    int SNP_count = 0 ;
    int solidSNP = 0;
    int InDel_count = 0 ;
    int SV_count = 0 ;
    std::ifstream ifF1(F1_VCF);
    if(!ifF1.is_open()) {
        std::cerr<<"open file "<<F1_VCF<<" failed !! exit ..."<<std::endl;
        exit(-1) ;
    }
    while(!std::getline(ifF1,line).eof()){
        if(!isValidVCF(line)) continue ;
        VI  temp;
        temp.InitFromVCF(line);
        variant_count ++ ;
        if(temp.type == V_Type::isSNP) SNP_count++;
        if(temp.type == V_Type::isInDel) InDel_count++;
        if(temp.type == V_Type::isSV) SV_count++;
        if( temp.type != V_Type::isSNP ) continue ;
        if( temp.htype == V_HType::type_1_1 ) continue ;
        if( temp.gt_str == "0/1" || temp.gt_str == "1/0" ||  temp.gt_str == "2/1" || temp.gt_str == "1/2" ) continue ;
        std::cout<<temp.ref_name<<'\t'<<temp.pos<<'\t'<<temp.alt1<<'\t'<<temp.alt2<<'\t'<<temp.phased_id<<'\n';
        solidSNP ++;
    }
    std::cerr<<"Loaded total\t"<<variant_count<<" variants from"<<F1_VCF<<std::endl;
    std::cerr<<"       SNPs\t"<<SNP_count<<std::endl;
    std::cerr<<"phased SNPs\t"<<solidSNP<<std::endl;
    std::cerr<<"       InDels\t"<<InDel_count<<std::endl;
    std::cerr<<"       SVs\t"<<SV_count<<std::endl<<std::endl;
    std::cerr<<"All done"<<std::endl;
}
