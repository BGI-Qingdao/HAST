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


namespace Test{
    int x ;
}

int main(int argc , char ** argv)
{
    // -------------------------------------------------------------

    // -------------------------------------------------------------
    std::BGIQD::APP::FileNames fNames;
    Test::x = 1 ;
    fNames.Init(prefix.to_string());

}
