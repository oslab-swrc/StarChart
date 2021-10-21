CPPFLAGS=-Ofast -march=native -flto -std=c++14 -fopenmp -fpic -Wpedantic -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare

ssa: src/code.c
	g++ $(CPPFLAGS) $(INCLUDEFLAGS) src/code.c -g -o  bin/ssa

clean: 
	rm -f bin/ssa
