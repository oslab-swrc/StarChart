CPPFLAGS=-Ofast -march=native -flto -std=c++14 -fopenmp -fpic -Wpedantic -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
INCLUDEFLAGS=-I src/lib/include -I src/lib/external/eigen -I src/lib/external/simple-serializer

ssa: src/code.cc
	g++ $(CPPFLAGS) $(INCLUDEFLAGS) src/code.cc -g -o  bin/ssa

clean: 
	rm -f bin/ssa
