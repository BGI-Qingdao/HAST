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
 Brief   : Merge two haplotype vcf and print all SNP in two alleles style.\n\
\n\
 Notice  : I only use chr1-chr22 and I will change CHROM from chr1 to 1.\n\
\n\
 Usage   : \n\
   ./MergeHapSNP  H1.vcf H2.vcf >paired_snp.txt\n\
\n\
 Output  :\n\
   CHROM\tPOS\tN1\tN2\n\
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
    bool isChr1_22() const {
        if( ref_name.size() >3  && ref_name.substr(0,3) == "chr" ){
            if( ref_name.size() <=5 && std::isdigit(ref_name.at(3)))
                return true ;
            else 
                return false;
        }
        else {
            return false ;
        }
    }
    std::string refName_no_chr()const {
        return ref_name.substr(3);
    }
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

void LoadChr1_22_SNP_FromVCF(const std::string & filename,  std::map<std::string , std::map<int , VI> > & amap)
{
    std::ifstream ifP1(filename);
    if(!ifP1.is_open()) {
        std::cerr<<"open file "<<filename<<" failed !! exit ..."<<std::endl;
        exit(-1) ;
    }
    std::string line;
    while(!std::getline(ifP1,line).eof()){
        if(!isValidVCF(line)) continue ;
        VI  temp;
        temp.InitFromVCF(line);
        if( temp.isChr1_22() && temp.type == V_Type::isSNP ) {
            assert(temp.gt_str == "1/1" ) ;
            amap[temp.ref_name][temp.pos] = temp ;
        }
    }
}

std::map<std::string , std::map<int , VI>  >  H1;
std::map<std::string , std::map<int , VI>  >  H2;

struct HapSNP {
    std::string ref_name ;
    int pos ;
    std::string alt1;
    std::string alt2;
};

std::map<std::string, std::map<int , HapSNP> > all_hapsnps;

bool isNewHAPSNP(const VI & vi) {
    if( all_hapsnps.find(vi.ref_name) == all_hapsnps.end()) return true ;
	const auto & chrs =  all_hapsnps.at(vi.ref_name);
    if( chrs.find(vi.pos) ==  chrs.end()) return true ;
    return false ;
}
void AddNewHapSNP( const VI & vi , std::string N2) {
    all_hapsnps[vi.ref_name][vi.pos].ref_name = vi.refName_no_chr();
    all_hapsnps[vi.ref_name][vi.pos].pos = vi.pos ;
    all_hapsnps[vi.ref_name][vi.pos].alt1 = vi.alt1 ;
    all_hapsnps[vi.ref_name][vi.pos].alt2 = N2 ;
}

void UpdateHapSNP( const VI & vi , std::map<std::string , std::map<int , VI>  >  another){
    if( ! isNewHAPSNP(vi) ) return ; 
    if( another.find(vi.ref_name) == another.end() ) {
        AddNewHapSNP(vi, vi.ref);
        return ;
    }
    const auto & chrs =  another.at(vi.ref_name);
    if( chrs.find(vi.pos) == chrs.end() ) {
        AddNewHapSNP(vi, vi.ref);
        return ;
    }
    const auto avi = chrs.at(vi.pos);
    AddNewHapSNP(vi,avi.alt1);
}

int main(int argc , char ** argv)
{
    if(argc != 3 ) { 
        PrintUsage();
        return 0;
    }
    std::string H1_VCF(argv[1]);
    std::string H2_VCF(argv[2]);
    if( H1_VCF == "-h" || H1_VCF == "--help" ){
        PrintUsage();
        return 0;
    }
    LoadChr1_22_SNP_FromVCF(H1_VCF,H1);
    LoadChr1_22_SNP_FromVCF(H2_VCF,H2);

    // Check all SNP in H1
    for(const auto & chroms_pair : H1) {
        const auto & chroms = chroms_pair.second;
        for( const auto & vi_pair : chroms ) {
            const auto & vi = vi_pair.second;
            UpdateHapSNP(vi,H2);
        }
    }
    // Check all SNP in H2
    for(const auto & chroms_pair : H2) {
        const auto & chroms = chroms_pair.second;
        for( const auto & vi_pair : chroms ) {
            const auto & vi = vi_pair.second;
            UpdateHapSNP(vi,H1);
        }
    }
    // print final result
    for(const auto & chroms_pair : all_hapsnps) {
        const auto & chroms = chroms_pair.second;
        for( const auto & vi_pair : chroms ) {
            const auto & vi = vi_pair.second;
            std::cout<<vi.ref_name<<'\t';
            std::cout<<vi.pos<<'\t';
            std::cout<<vi.alt1<<'\t';
            std::cout<<vi.alt2<<'\n';
        }
    }
}
