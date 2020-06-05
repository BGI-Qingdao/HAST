
.PHONY:all clean

all:
	cd 01.classify_stlfr_reads && make && cd ../
	cd 03.mkoutput_by_fabulous2.0/src_main && make && cd ../../

clean:
	cd 01.classify_stlfr_reads && make clean && cd ../
	cd 03.mkoutput_by_fabulous2.0/src_main && make clean && cd ../../
