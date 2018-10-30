all: districting

districting: *.cpp *.h
	g++ -Wall graph.cpp main.cpp io.cpp -o districting -I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -lgurobi81 -lgurobi_g++4.2 -std=c++11 -O2 -Wno-sign-compare -g
