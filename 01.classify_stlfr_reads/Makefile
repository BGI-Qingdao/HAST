.PHONY: all tool clean

all: classify tool 

classify : classify.cpp gzstream/gzstream.C gzstream/gzstream.h kmer/kmer.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11  classify.cpp gzstream.o -lz -lpthread -o classify
	#g++ -g -std=c++11 -static classify.cpp gzstream.o -lz -lpthread -o classify

clean :
	rm classify *.o gzstream/*.a
