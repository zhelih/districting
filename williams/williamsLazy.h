#pragma once
#ifndef LAZYCONSTRAINTS_H
#define LAZYCONSTRAINTS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "GRBInterface.h"
#include "gurobi_c++.h"
#include <ctime>
#include <memory>
#include <sstream>
#include "KGraph.h"

extern double copy_time;
using namespace std;

class LazyConstraints : public GRBCallback //Lazy constraints for redistricting problem
{
private:
	long** x;
	long *r;

public:
	LazyConstraints(GRBVar **xvars, GRBVar *rvars, KGraph &g, KGraph &d, long &K, long &L, long &U);
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;

};
#endif
