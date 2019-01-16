all: districting

districting: *.cpp *.h version.c
	g++ -Wall -Wextra version.c LazyConstraints.cpp graph.cpp main.cpp io.cpp hess.cpp flow.cpp cut.cpp diff_ul.cpp -o districting -I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -lgurobi81 -lgurobi_g++4.2 -std=c++11 -O2 -Wno-sign-compare -g

version.c: ../.git/HEAD ../.git/index
	echo "const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

gridgen: gridgen.cpp
	g++ -Wall -Wextra graph.cpp gridgen.cpp -o gridgen -std=c++11

clean:
	rm ./districting
	rm ./version.c
