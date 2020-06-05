#ifndef __APPCOMMON_IDX_H__
#define __APPCOMMON_IDX_H__ 

#include <vector>
#include <map>
#include <sstream>
namespace BGIQD {
    namespace APP {
        struct Idx {
            unsigned int scaffold_id ;
            std::vector<int> indexs;
            bool is_valid() const { 
                return indexs.size() >1 && indexs.size()%2==0 ;
            }
            bool is_single() const {
                return indexs.size() == 2 ;
            }
            bool is_multi() const {
                return indexs.size() >2 ;
            }
            std::map<int,int> phase_parts() const { 
                std::map<int,int> ret ;
                if( !is_valid() || is_single() ) return ret ;
                for( int i = 1 ; i < (int)indexs.size()-2 ; i+= 2 ) {
                    ret[indexs.at(i)] = indexs.at(i+1);
                }
                return ret ;
            }
            std::map<int,int> homo_parts() const { 
                std::map<int,int> ret ;
                if( !is_valid() ) return ret ;
                for( int i = 0 ; i <= (int)indexs.size()-2 ; i+= 2 ) {
                    ret[indexs.at(i)] = indexs.at(i+1);
                }
                return ret ;
            }
            void InitFromString(const std::string & line ) {
                std::istringstream ist(line);
                ist>>scaffold_id;
                while(!ist.eof()){
                    int index ;
                    ist >> index ;
                    indexs.push_back(index) ;
                }
            }
            std::string ToString() const {
                std::ostringstream ost ;
                ost<<scaffold_id;
                for(auto x : indexs ) ost<<' '<<x;
                return ost.str() ;
            }
        };
    }
}
#endif //__APPCOMMON_IDX_H__
