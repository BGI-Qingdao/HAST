#ifndef __APPCOMMON_SEGMENTFA_H__
#define __APPCOMMON_SEGMENTFA_H__

#include <string>
#include <sstream>
#include <cassert>

namespace BGIQD { 
    namespace APP {
        struct Scaff_Seg_Head
        {
            unsigned int scaff_id;

            int  seq_index ; // index start from 0 
                             // 0 : 0-1 seg 
                             // 1 : 1-2 seg 
                             //-1 : invalid

            int  phase_id ;  // 0 : homo sequence
                             // 1 : supernova phase group 1
                             // 2 : supernova phase group 2
                             //-1 : invalid

            void Init( const std::string & line ) {
                assert( ! line.empty() ) ;
                if( line[0] == '>' )
                    sscanf(line.c_str() , ">%u_%d_%d",&scaff_id,&seq_index,&phase_id);
                else 
                    sscanf(line.c_str() , "%u_%d_%d",&scaff_id,&seq_index,&phase_id);
            }
            bool Valid() const { 
                return seq_index != -1 && phase_id != -1  && scaff_id != (unsigned int)-1 ;
            }
            std::string Head() const { 
                assert( seq_index != -1 && phase_id != -1  && scaff_id != -1 );
                return ">" 
                    + std::to_string(scaff_id)
                    + "_" 
                    + std::to_string(seq_index) 
                    + "_"
                    + std::to_string(phase_id);
            }

            void Reset() { scaff_id = -1 ; phase_id = -1 ; seq_index = -1; } 
        };
    }
}

#endif
