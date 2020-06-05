#ifndef __APPCOMMON_FILENAME_H__
#define __APPCOMMON_FILENAME_H__

#include <string>
#include "common/string/stringtools.h"

namespace BGIQD {
    namespace APP {

        struct FileNames
        {
            void Init(const std::string & prefix)
            {
                m_prefix = prefix ;
            }

#define DEFINE_SUFFIX(name,suffix) \
            std::string name() const { return m_prefix + suffix ; } \
            std::string name(int round) const \
            {\
                if ( round == 0 )\
                {\
                    return name();\
                }\
                else\
                {\
                    return m_prefix + suffix +"_round_"+BGIQD::STRING::itos(round) ;\
                }\
            }\
            std::string name(std::string middle) const \
            {\
                if( middle == "" )\
                {\
                    return name() ;\
                }\
                else\
                {\
                    return m_prefix + "." + middle+  suffix ;\
                }\
            }\

            //Split output 
            DEFINE_SUFFIX(phb1_fa,".phb.1.fa");
            DEFINE_SUFFIX(phb2_fa,".phb.2.fa");
            DEFINE_SUFFIX(homo_fa,".homo.fa");
            //Merge input
            DEFINE_SUFFIX(homo_ids,".homo.ids");
            DEFINE_SUFFIX(father_ids,".father.ids");
            DEFINE_SUFFIX(mother_ids,".mother.ids");
            //Merge output
            DEFINE_SUFFIX(merge_father_ids,".merge.father.ids");
            DEFINE_SUFFIX(merge_mother_ids,".merge.mother.ids");
            //Gen output
            DEFINE_SUFFIX(father_fa,".father.fa");
            DEFINE_SUFFIX(mother_fa,".mother.fa");
            DEFINE_SUFFIX(father_idx,".father.idx");
            DEFINE_SUFFIX(mother_idx,".mother.idx");
            private:
                std::string m_prefix;
        };
    }//namespace APP
}//namespace BGIQD

#endif //__SOAP2_FILENAME_H__
