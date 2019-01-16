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


class LazyConstraints : public GRBCallback //Lazy constraints for solving CUT1 in political redistricting
{
private:
  GRBVar **varsx;
  GRBVar **varsy;
  graph* g1;
public:
	LazyConstraints(GRBVar **xvars, GRBVar **yvars, graph *g);
  virtual ~LazyConstraints() {}
	long numCallbacks;
	double TotalCallbackTime;
	long numLazyCuts;
protected:
	void callback();
};

#endif
