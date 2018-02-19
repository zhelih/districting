all:
	g++ -Wall graph.cpp main.cpp -o gerry_assign -I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -lgurobi70 -lgurobi_g++4.1 -laes70 -std=c++11 -O2
