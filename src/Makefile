INTELROOT=/home/lykhovyd/intel
MKLROOT=$(INTELROOT)/mkl
MKLLIB=$(MKLROOT)/lib/intel64_lin
MKLINCLUDE=$(MKLROOT)/include


MKL_FLAGS=-L$(MKLLIB) -I$(MKLINCLUDE)/intel64/lp64 -lmkl_blas95_lp64 -Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
GUROBI_FLAGS=-I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -L"$(INTELROOT)/lib/intel64"  -lgurobi81 -lgurobi_g++5.2
GENERAL_FLAGS=-std=c++11 -O3 -Wno-sign-compare -Wall -Wextra


# to save some space
COMMON_OBJ=version.c graph.o lagrange.o io.o hess.o flow.o cut.o ralg.o


all: districting ralg_hot_start translate

districting: $(COMMON_OBJ) main.o
	g++ main.o $(COMMON_OBJ) -o districting $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKL_FLAGS)
	cp districting ../

ralg_hot_start: version.c graph.o lagrange.o main_hot_start.o io.o hess.o flow.o cut.o ralg.o
	g++ main_hot_start.o $(COMMON_OBJ) -o ralg_hot_start $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKL_FLAGS)
	cp ralg_hot_start ../

ralg.o: ralg/*
	g++ -c -Wall ralg/ralg.cpp -O2 -I"$(MKLINCLUDE)" -std=c++11

%.o: %.cpp *.h
	g++ $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKL_FLAGS) -c $< -o $@

version.c: ../.git/HEAD ../.git/index
	echo "const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

gridgen: gridgen.cpp
	g++ -Wall -Wextra graph.cpp gridgen.cpp -o gridgen -std=c++11

translate: translate.cpp
	g++ -Wall -Wextra translate.cpp -o translate -std=c++11

clean:
	rm ./districting
	rm ./ralg_hot_start
	rm ./translate
	rm *.o
	rm ./version.c
