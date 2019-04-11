all: districting

districting: *.cpp *.h version.c
	g++ -Wall -Wextra version.c graph.cpp lagrange.cpp main.cpp io.cpp hess.cpp flow.cpp cut.cpp diff_ul.cpp -o districting -I"$(GUROBI_HOME)/include" -L"$(GUROBI_HOME)/lib" -lgurobi_c++ -lgurobi81 -std=c++11 -O2 -Wno-sign-compare -g

version.c: ../.git/HEAD ../.git/index
	echo "const char *gitversion = \"$(shell git describe --tags --always)\";" > $@

gridgen: gridgen.cpp
	g++ -Wall -Wextra graph.cpp gridgen.cpp -o gridgen -std=c++11

translate: translate.cpp
	g++ -Wall -Wextra translate.cpp -o translate -std=c++11

clean:
	rm ./districting
	rm ./version.c
