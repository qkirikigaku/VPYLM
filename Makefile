#!bin/make

CXXFLAGS = -O2 -Wall

compile:
	g++ -O2 -Wall src/main.cpp -o VPYLM -std=c++11

debug:
	g++ -O0 -Wall -g src/main.cpp -o VPYLM_debug -std=c++11

profile:
	g++ -g -O2 src/main.cpp -o VPYLM_profile -std=c++11
	iprofiler -timeprofiler -T 300 ./VPYLM_profile test Upstream 1

test:
	make compile
	./VPYLM test

test_debug:
	make debug
	lldb VPYLM_debug test
