

classify : classify.cpp
	g++ -c  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -std=c++11 classify.cpp gzstream.o -lz -lpthread -o classify
