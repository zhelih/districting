#include <vector>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k);

