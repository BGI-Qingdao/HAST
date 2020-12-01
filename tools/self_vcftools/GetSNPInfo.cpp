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
		Brief   : Get only SNP information from a VCF file.\n\
		\n\
		Usage   : \n\
		./GetSNPInfo XXX.vcf \n\
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

enum V_Type {
	isSNP = 0,
	isInDel = 1 ,
	isSV = 2
};

// Variant Information
struct VI {
	std::string ref_name;
	int pos ;
	std::string ref;
	std::string alt;
	std::string alt2;
	std::string alt1;
	std::set<std::string> seqs;
	std::string  gt_str;
	V_HType htype;
	V_Type type ;
	// Func1
	void InitFromVCF(const std::string & line){
		auto items = split(line,'\t');
		ref_name=items[0];
		pos = std::stoi(items[1]);
		ref=items[3];
		alt=items[4];
		auto v_alts=split(items[4],',');
		for( auto x : v_alts ){
			seqs.insert(x);
		}
		int GT_index=-1; int PS_index=-1;
		auto describe = split(items[8],':');
		for( int i = 0 ; i < describe.size() ; i++ ){
			if( describe[i] == "GT" ) GT_index=i;
			if( describe[i] == "PS" ) PS_index=i;
		}
		auto datas=split(items[9],':');
		if( GT_index >= 0 )
		{
		gt_str = datas[GT_index];
		if( gt_str  == "0|1" || gt_str== "0/1" ||  gt_str  == "1|0" ||  gt_str == "1/0" ) {
			htype =V_HType::type_0_1;
		}
		else if( gt_str == "1|1" || gt_str == "1/1" ) {
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
	if( std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help" ) {
		PrintUsage();
		return 0;
	}
	std::string A_VCF(argv[1]);
	std::ifstream ifA(A_VCF);
	if(!ifA.is_open()) {
		std::cerr<<"open file failed !! exit ..."<<std::endl;
		return 0 ;
	}
	// Load A
	std::string line;
	int variant_count = 0;
	int SNP_count = 0 ;
	int InDel_count = 0 ;
	int SV_count = 0 ;
	while(!std::getline(ifA,line).eof()){
		if(!isValidVCF(line)) continue ;
		VI  temp;
		temp.InitFromVCF(line);
		variant_count ++ ;
		if(temp.type == V_Type::isInDel) InDel_count++;
		if(temp.type == V_Type::isSV) SV_count++;
		if(temp.type == V_Type::isSNP){
			SNP_count++;
			std::cout<<temp.ref_name<<'\t';
			std::cout<<temp.pos<<'\t';
			std::cout<<temp.alt1<<'\t';
			std::cout<<temp.alt2<<'\n';
		}
	}
	std::cerr<<"Loaded total\t"<<variant_count<<" variants from"<<argv[1]<<std::endl;
	std::cerr<<"       SNPs\t"<<SNP_count<<std::endl;
	std::cerr<<"       InDels\t"<<InDel_count<<std::endl;
	std::cerr<<"       SVs\t"<<SV_count<<std::endl;
	std::cerr<<"All done"<<std::endl;
	return 0;
}
