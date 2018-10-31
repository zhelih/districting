#ifndef _MODELS_H
#define _MODELS_H

#include <vector>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

// build hess model and return x variables
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k);
// add SCF constraints to model with hess variables x
void build_scf(GRBModel* model, GRBVar** x, graph* g);
// add MCF2 constraints to model with hess variables x
void build_mcf(GRBModel* model, GRBVar** x, graph* g);
// add CUT2 constraints to model with hess variables x
void build_cut2(GRBModel* model, GRBVar** x, graph* g);
#endif
