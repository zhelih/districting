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
// add MCF constraints to model with hess variables x
void build_mcf1(GRBModel* model, GRBVar** x, graph* g);
void build_mcf2(GRBModel* model, GRBVar** x, graph* g);
// add CUT constraints to model with hess variables x (lazy)
void build_cut1(GRBModel* model, GRBVar** x, graph* g);
void build_cut2(GRBModel* model, GRBVar** x, graph* g);

// add UL1 & UL2 instance and return x variables
GRBVar** build_UL_1(GRBModel* model, graph* g, const vector<int>& population, int k);
GRBVar** build_UL_2(GRBModel* model, graph* g, const vector<int>& population, int k);

#endif
