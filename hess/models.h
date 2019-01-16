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

class Cut2Callback: public GRBCallback
{
  // memory for a callback
  private:
    double **x; // x values
    int* visited; // dfs marks
    int* aci; // A(C_i) set
    GRBVar** grb_x; // x variables
    graph* g;
    int n;
    std::vector<int> s; // stack for DFS
  public:
    int numCallbacks; // number of callback calls
    double callbackTime; // cumulative time in callbacks

    Cut2Callback(GRBVar** grb_x_, graph *g_) : grb_x(grb_x_), g(g_), numCallbacks(0), callbackTime(0.)
    {
      n = g->nr_nodes;
      x = new double*[n];
      for(int i = 0; i < n; ++i)
        x[i] = new double[n];
      visited = new int[n];
      aci = new int[n];
      s.reserve(n);
    }
    virtual ~Cut2Callback()
    {
      delete [] aci;
      delete [] visited;
      for(int i = 0; i < n; ++i)
        delete [] x[i];
      delete [] x;
    }
  protected:
    void callback();
};

Cut2Callback* build_cut2(GRBModel* model, GRBVar** x, graph* g);

// add UL1 & UL2 instance and return x variables
GRBVar** build_UL_1(GRBModel* model, graph* g, const vector<int>& population, int k);
GRBVar** build_UL_2(GRBModel* model, graph* g, const vector<int>& population, int k);

#endif
