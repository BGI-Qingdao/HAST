.PHONY: all

all: HapInherit	GetSNPOnly 	InheritSolidSNP	PhasedSNP MergeHapSNP CalcHD

HapInherit: HapInherit.cpp
	g++ -g -Werror -std=c++11 HapInherit.cpp -o HapInherit

GetSNPOnly:	GetSNPOnly.cpp
	g++ -g -Werror -std=c++11 GetSNPOnly.cpp -o GetSNPOnly

InheritSolidSNP:	InheritSolidSNP.cpp
	g++ -g -Werror -std=c++11 InheritSolidSNP.cpp -o InheritSolidSNP

PhasedSNP: PhasedSNP.cpp
	g++ -g -Werror -std=c++11 PhasedSNP.cpp -o PhasedSNP

MergeHapSNP : MergeHapSNP.cpp
	g++ -g -Werror -std=c++11 MergeHapSNP.cpp -o MergeHapSNP

CalcHD : CalcHD.cpp
	 g++ -g -Werror -std=c++11 CalcHD.cpp -o CalcHD
