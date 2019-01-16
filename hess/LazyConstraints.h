#ifndef LAZYCONSTRAINTS_H
#define LAZYCONSTRAINTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include "models.h"
#include "gurobi_c++.h"
#include <ctime>
#include <memory>
#include <sstream>
#include "graph.h"


extern double copy_time;
using namespace std;

class LazyConstraints : public GRBCallback //Lazy constraints for solving CUT1 in political redistricting
{
public:
	LazyConstraints(GRBVar **xvars, GRBVar **yvars, graph *g);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;
};

#endif
