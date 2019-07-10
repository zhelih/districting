MKLROOT=/home/lykhovyd/intel/mkl
MKLLIB=$(MKLROOT)/lib/intel64_lin
MKLINCLUDE=$(MKLROOT)/include

all: districting

districting: *.cpp *.h version.c ralg.o
	g++ -Wall -Wextra version.c graph.cpp lagrange.cpp main.cpp io.cpp hess.cpp flow.cpp cut.cpp preprocess.cpp -o districting -I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib"  -lgurobi81 -lgurobi_g++5.2 -std=c++11 -O2 -Wno-sign-compare ralg.o -L$(MKLLIB) -I$(MKLINCLUDE)/intel64/lp64 -lmkl_blas95_lp64 -Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

version.c: ../.git/HEAD ../.git/index
	echo "const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

gridgen: gridgen.cpp
	g++ -Wall -Wextra graph.cpp gridgen.cpp -o gridgen -std=c++11

translate: translate.cpp
	g++ -Wall -Wextra translate.cpp -o translate -std=c++11

clean:
	rm ./districting
	rm ./version.c

ralg.o: ralg/*
	g++ -c -Wall ralg/ralg.cpp -O2 -I"$(MKLINCLUDE)" -std=c++11
