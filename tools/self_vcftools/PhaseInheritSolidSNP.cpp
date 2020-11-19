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
 Brief   : Find the solid SNP inherited from two diploid parents to a diploid offspring.\n\
           solid refer to the two nucleotides of SNP can be reconstructed from parents\n\
           without any mutation.\n\
\n\
 Usage   : \n\
   ./HapInheritSolidSNP HG003.vcf HG004.vcf HG002.vcf  >solid_snp.vcf \n\
\n\
 Output  :\n\
   for each SNP in offspring VCF , it will print :\n\
   CHROM\tPOS\t.\tREF\tALT\tGT\tSGT	PIT	MIT\n\
\n\
 Detail  :\n\
   the columns one to five are the same as offspring VCF.\n\
   SGT refer to short-GT:\n\
       * 0_1 represent { 0/1 , 0|1 , 1|0 , 1/0 }\n\
       * 1_1 represent { 1/1 , 1|1 } \n\
       * 1_2 represent { 1/2 , 1|2 , 2|1 , 2/1 }\n\
\n\
   (P/M)IT refer to inherited type of (p/m)aternal:\n\
       * 0 represent variant not happened in parent.\n\
       * 1 represent parent has variant at the same position\n\
           but inherit reference allele to offspring.\n\
       * 2 represent parent has variant at the same position\n\
           and inherit the variant allele to offspring as offspring alt1.\n\
       * 3 represent parent has variant at the same position\n\
           and inherit the variant allele to offspring as offspring alt2.\n\
       * 4 represent parent has the same alleles as offsprint\n\
           therefor I can not make decision without another\n\
           parental inforamtion.\n\
       * 5 represent parent has variant at the same position\n\
           but none of the alleles is the same as any of the\n\
           offspring alleles. \n\
           (Indel and SV may hit this because I use exact match)\n";
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


std::map<std::string , std::map<int , VI> > P1_map;
std::map<std::string , std::map<int , VI> > P2_map;
std::map<std::string , std::map<int , VI> > F1_map;

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

bool isSNPSolid(const VI & vi , A_in_B_Type p1 , A_in_B_Type p2)
{
    if( vi.type != V_Type::isSNP ) return false ;
    if( vi.htype == V_HType::type_1_1 ) {
       if(p1 == A_in_B_Type::A_in_B_alt1 && p2 == A_in_B_Type::A_in_B_alt1 ) 
           return true;
       else 
           return false;
    }
    else if ( vi.htype ==  V_HType::type_0_1 ) {
        if( p1 == A_in_B_Type::A_not_in_B || p1 == A_in_B_Type::A_in_B_ref ) {
            if( p2 ==  A_in_B_Type::A_in_B_alt1 || p2  == A_in_B_Type::A_in_B_all ) 
                return true ;
            else 
                return false;
        }
        else if ( p1 == A_in_B_Type::A_in_B_alt1 ) {
            if( p2 == A_in_B_Type::A_in_B_all 
                    ||p2 == A_in_B_Type::A_not_in_B
                    ||p2 ==  A_in_B_Type::A_in_B_ref )
                return true ;
            else 
                return false;
        }
        else if  ( p1 == A_in_B_Type::A_in_B_all ) {
            if( p2 == A_in_B_Type::A_in_B_all 
                    ||p2 == A_in_B_Type::A_not_in_B
                    ||p2 == A_in_B_Type::A_in_B_alt1
                    ||p2 ==  A_in_B_Type::A_in_B_ref )
                return true ;
            else 
                return false;
        }
        else 
            return false ;
    } else if ( vi.htype == V_HType::type_1_2 ) {
        if( p1 == A_in_B_Type::A_in_B_all ){
            if( p2 == A_in_B_Type::A_in_B_all 
                || p2 == A_in_B_Type::A_in_B_alt1 
                || p2 == A_in_B_Type::A_in_B_alt2 )
                return true ;
            else
                return false ;
        }
        else if ( p1 == A_in_B_Type::A_in_B_alt1 ) {
            if( p2 == A_in_B_Type::A_in_B_all 
                || p2 == A_in_B_Type::A_in_B_alt2 )
                return true ;
            else
                return false ;
        }
        else if ( p1 == A_in_B_Type::A_in_B_alt2 ) {
            if( p2 == A_in_B_Type::A_in_B_all 
                || p2 == A_in_B_Type::A_in_B_alt1 )
                return true ;
            else
                return false ;
        }
        else
            return false;
    }
    return false;
}

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
    LoadFromVCF(F1_VCF,F1_map);
    // Check all SNP in F1
    for(const auto & chroms_pair : F1_map) {
        const auto & chroms = chroms_pair.second;
        for( const auto & vi_pair : chroms ) {
            const auto & vi = vi_pair.second;
            if( vi.type != V_Type::isSNP ) continue ;
            A_in_B P1 ;
            A_in_B P2 ;
            P1.type = V_in_Parent(vi,P1_map,P1.inherit);
            P2.type = V_in_Parent(vi,P2_map,P2.inherit);
        }
    }
    std::cerr<<"All done"<<std::endl;
}
