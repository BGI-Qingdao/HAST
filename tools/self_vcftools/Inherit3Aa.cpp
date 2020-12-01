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
 Brief   : Find the 3Aa SNP inherited from two diploid parents to a diploid offspring.\n\
\n\
 Usage   : \n\
   ./Inherit3Aa HG003.vcf HG004.vcf HG002.vcf  >solid_snp.vcf \n\
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
enum A_in_B_Type{
    A_not_in_B = 0 ,
    A_in_B_ref = 1 ,
    A_in_B_alt1 = 2 ,
    A_in_B_alt2 = 3 ,
    A_in_B_all = 4 , 
    A_diff_B = 5 ,
};

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
        if( datas[0] == "0|1" || datas[0] == "0/1" ||  datas[0] == "1|0" ||  datas[0] == "1/0" ) {
            htype =V_HType::type_0_1;
        }
        else if( datas[0] == "1|1" || datas[0] == "1/1" ) {
            htype =V_HType::type_1_1;
        } else {
            htype =V_HType::type_1_2;
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

struct A_in_B {
    A_in_B_Type type ;
    std::string inherit;
    //CHROM	POS	REF	ALT	GT	SGT	VT	IT	IS
    void Print(const VI & vi) const {
        std::cout<<vi.ref_name<<'\t'<<vi.pos<<'\t'<<vi.ref<<'\t'<<vi.alt<<'\t';
        std::cout<<vi.gt_str<<'\t';
        std::cout<<V_HType_ToStr(vi.htype)<<'\t';
        std::cout<<V_Type_ToStr(vi.type)<<'\t';
        std::cout<<int(type)<<'\t';
        if( type == A_in_B_Type::A_not_in_B ) {
            std::cout<<"."<<'\n';
        } else if ( type == A_in_B_Type::A_in_B_ref 
                || type == A_in_B_Type::A_in_B_alt1 
                || type == A_in_B_Type::A_in_B_alt2 ) {
            std::cout<<inherit<<'\n';
        } else if ( type == A_in_B_Type::A_in_B_all ) {
            std::cout<<"*"<<'\n';
        } else if ( type ==  A_in_B_Type::A_diff_B ) {
            std::cout<<"."<<'\n';
        }
    }
};

bool isValidVCF(const std::string & line) {
    return line.size()>0 && line.at(0) != '#' ;
}

A_in_B_Type V_in_Parent( const VI & item , const std::map<std::string , std::map<int , VI> > & parent , std::string & c ) {
    if( parent.find(item.ref_name) == parent.end() ) return A_in_B_Type::A_not_in_B ;
    const auto & chrs =  parent.at(item.ref_name);
    if( chrs.find(item.pos) == chrs.end() ) return A_in_B_Type::A_not_in_B ;
    const auto & vi = chrs.at(item.pos) ;
    int match_num = 0 ;
    std::string match_str;
    int i=0;
    for(const auto & seq : item.seqs ) {
        i++;
        if( vi.seqs.find(seq) != vi.seqs.end() ){
            match_num +=i ;
            match_str = seq;
        }
    }
    if( match_num == 0 ) return A_in_B_Type::A_diff_B;
    if( match_num == 3 ) return A_in_B_Type::A_in_B_all ;
    if( match_num >3 ) { assert(0) ; return A_in_B_Type::A_in_B_all; }
    if( match_num == 1 || match_num == 2) {
        c = match_str;
        if( item.htype == V_HType::type_0_1 ) {
            if( match_str == item.ref )
                return A_in_B_Type::A_in_B_ref ;
            else 
                return A_in_B_Type::A_in_B_alt1 ;
        } 
        else if ( item.htype == V_HType::type_1_1 ) {
            return A_in_B_Type::A_in_B_alt1;
        }
        else{
            if( match_num == 1 )
                return A_in_B_Type::A_in_B_alt1;
            else 
                return A_in_B_Type::A_in_B_alt2;
        }
    }else { assert(0) ; return  A_in_B_Type::A_in_B_ref ; }

}



void LoadFromVCF(const std::string & filename,  std::map<std::string , std::map<int , VI> > & amap)
{
    std::ifstream ifP1(filename);
    if(!ifP1.is_open()) {
        std::cerr<<"open file "<<filename<<" failed !! exit ..."<<std::endl;
        exit(-1) ;
    }
    std::string line;
    int variant_count = 0;
    int SNP_count = 0 ;
    int InDel_count = 0 ;
    int SV_count = 0 ;
    while(!std::getline(ifP1,line).eof()){
        if(!isValidVCF(line)) continue ;
        VI  temp;
        temp.InitFromVCF(line);
        variant_count ++ ;
        if(temp.type == V_Type::isSNP) SNP_count++;
        if(temp.type == V_Type::isInDel) InDel_count++;
        if(temp.type == V_Type::isSV) SV_count++;
        amap[temp.ref_name][temp.pos] = temp ;
    }
    std::cerr<<"Loaded total\t"<<variant_count<<" variants from"<<filename<<std::endl;
    std::cerr<<"       SNPs\t"<<SNP_count<<std::endl;
    std::cerr<<"       InDels\t"<<InDel_count<<std::endl;
    std::cerr<<"       SVs\t"<<SV_count<<std::endl<<std::endl;
}

bool isSNP3Aa(const VI & vi , A_in_B_Type p1 , A_in_B_Type p2)
{
    if( vi.type != V_Type::isSNP ) return false ;
    if( vi.htype == V_HType::type_0_1 ) {
       if(p1 == A_in_B_Type::A_in_B_all && p2 == A_in_B_Type::A_in_B_all ) 
           return true;
       else 
           return false;
    }
    return false ;
}

std::map<std::string , std::map<int , VI> > P1_map;
std::map<std::string , std::map<int , VI> > P2_map;
std::map<std::string , std::map<int , VI> > F1_map;
int main(int argc , char ** argv)
{
    if(argc != 4 ) { 
        PrintUsage();
        return 0;
    }
    std::string P1_VCF(argv[1]);
    std::string P2_VCF(argv[2]);
    std::string F1_VCF(argv[3]);

    LoadFromVCF(P1_VCF,P1_map);
    LoadFromVCF(P2_VCF,P2_map);
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
        A_in_B P1 ;
        A_in_B P2 ;
        P1.type = V_in_Parent(temp,P1_map,P1.inherit);
        P2.type = V_in_Parent(temp,P2_map,P2.inherit);
        if( isSNP3Aa(temp,P1.type,P2.type) )
        {
            std::cout<<line<<std::endl;
            solidSNP ++ ;
        }
    }
    std::cerr<<"Loaded total\t"<<variant_count<<" variants from"<<F1_VCF<<std::endl;
    std::cerr<<"       SNPs\t"<<SNP_count<<std::endl;
    std::cerr<<" solid SNPs\t"<<solidSNP<<std::endl;
    std::cerr<<"       InDels\t"<<InDel_count<<std::endl;
    std::cerr<<"       SVs\t"<<SV_count<<std::endl<<std::endl;
    std::cerr<<"All done"<<std::endl;
}
