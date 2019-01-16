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
class HessCallback : public GRBCallback
{
protected:
  GRBVar** grb_x; // x variables
  double** x; // x values
  graph* g; // graph pointer
  int n; // g->nr_nodes
public:
  int numCallbacks; // number of callback calls
  double callbackTime; // cumulative time in callbacks
  int numLazyCuts;
  HessCallback(GRBVar** grb_x_, graph* g_) : grb_x(grb_x_), g(g_), numCallbacks(0), callbackTime(0.), numLazyCuts(0)
  {
    n = g->nr_nodes;
    x = new double*[n];
    for(int i = 0; i < n; ++i)
      x[i] = new double[n];
  }
  virtual ~HessCallback()
  {
    for(int i = 0; i < n; ++i)
      delete [] x[i];
    delete [] x;
  }
protected:
  void populate_x()
  {
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        x[i][j] = getSolution(grb_x[i][j]);
  }
};

// @return callback for delete only
HessCallback* build_cut1(GRBModel* model, GRBVar** x, graph* g);
HessCallback* build_cut2(GRBModel* model, GRBVar** x, graph* g);

// add UL1 & UL2 instance and return x variables
GRBVar** build_UL_1(GRBModel* model, graph* g, const vector<int>& population, int k);
GRBVar** build_UL_2(GRBModel* model, graph* g, const vector<int>& population, int k);

#endif
