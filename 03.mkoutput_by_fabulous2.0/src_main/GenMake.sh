#!/bin/bash

apps="\
 Split \
 MergePhaseResult \
 GenSq \
 classify \
"

jobs_o=" "

function GenApp()
{
local AppName=$1
jobs_o="$jobs_o \${$AppName""_o}"


echo """

$AppName"_cpp 	=	"$AppName".cpp"
$AppName"_o   =	"$AppName".o"
$AppName" : clean \${"$AppName"_o} \${source_o} ../bin"
	\${CXX} \${$AppName"_o} \${source_o} \${DEUBG_CXX}  -o "$AppName
	mv \$@ ../bin/

""">>Makefile

}

echo ".PHONY: all clean bin" >Makefile
echo """
CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11\\
				-I../\\
				-lz\\
				-lpthread\\

DEUBG_CXX  =	\${CXXFLAGS} -g
RELEASE_CXX=	\${CXXFLAGS}

source_cpp =	../common/files/file_reader.cpp \\
		   		../common/files/file_writer.cpp \\
		   		../common/files/gzstream.cpp \\
				../common/log/log.cpp\\
				../common/log/logfilter.cpp\\
				../common/time/timetools.cpp\\
				../common/string/stringtools.cpp\\
				../common/args/argsparser.cpp\\
				../biocommon/fasta/fasta.cpp\\
				../biocommon/seq/seq.cpp\\

source_o		= \${source_cpp:%.cpp=%.o}

.cpp.o:
	\${CXX} \${DEUBG_CXX} -c \$< -o \$@

jobs =$apps

all :  \${jobs}

""" >>Makefile
for x in $apps
do
    GenApp $x
done

echo "jobs_o=$jobs_o">>Makefile
echo """
dirty	   =\${jobs_o} \${jobs} \${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf \${dirty}
""" >>Makefile 
