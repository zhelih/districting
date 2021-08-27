GUROBI_FLAGS=-I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -lgurobi91 -lgurobi_g++5.2
GENERAL_FLAGS=-std=c++11 -O3 -Wno-sign-compare -Wall -Wextra


# to save some space
COMMON_OBJ=version.c graph.o lagrange.o io.o hess.o flow.o cut.o ralg.o
TARGETS=districting ralg_hot_start translate gridgen

all: check-env check-mkl-env $(TARGETS)

no-mkl: MKLINCLUDE=-Iralg/cblas
no-mkl: MKL_FLAGS= ralg/cblas/mkl_cblas.o
no-mkl: check-env ralg/cblas/mkl_cblas.o $(TARGETS)

districting: $(COMMON_OBJ) main.o
	g++ main.o $(COMMON_OBJ) -o districting $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKL_FLAGS)
	cp districting ../

ralg_hot_start: $(COMMON_OBJ) main_hot_start.o
	g++ main_hot_start.o $(COMMON_OBJ) -o ralg_hot_start $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKL_FLAGS)
	cp ralg_hot_start ../

ralg.o: ralg/*
	g++ -c -Wall ralg/ralg.cpp -O3 $(MKLINCLUDE) -std=c++11

%.o: %.cpp *.h
	g++ $(GENERAL_FLAGS) $(GUROBI_FLAGS) $(MKLINCLUDE) -c $< -o $@

update_version: ../.git/HEAD ../.git/index
	echo "const char *gitversion = \"$(shell git describe --tags --always)\";" > version.c

gridgen: gridgen.cpp graph.o
	g++ $(GENERAL_FLAGS) graph.o gridgen.cpp -o gridgen

translate: translate.cpp
	g++ $(GENERAL_FLAGS) translate.cpp -o translate

clean:
	rm -f ./districting
	rm -f ./ralg_hot_start
	rm -f ./translate
	rm -f *.o
	rm -f ./ralg/cblas/mkl_cblas.o

check-env:
ifndef GUROBI_HOME
  $(error GUROBI_HOME is undefined, please update your environment variables)
endif

check-mkl-env:
ifndef INTELROOT
	$(error INTELROOT is undefined, please use "make no-mkl" to use internal substitutes)
else
  MKLROOT:=$(INTELROOT)/mkl
  MKLLIB:=$(MKLROOT)/lib/intel64_lin
  MKLINCLUDE:=-I$(MKLROOT)/include -I$(MKLINCLUDE)/intel64/lp64
  MKL_FLAGS:=-L$(MKLLIB) -L$(INTELROOT)/lib/intel64 -lmkl_blas95_lp64 -Wl,--start-group $(MKLLIB)/libmkl_intel_lp64.a $(MKLLIB)/libmkl_intel_thread.a $(MKLLIB)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
endif

.PHONY: check-env check-mkl-env update_version
