#ifndef _MODELS_H
#define _MODELS_H

#include <vector>
#include "common.h"
#include "graph.h"
#include "io.h"
#include "gurobi_c++.h"

using namespace std;

typedef const vector<vector<bool>> cvv;

//auxilliary procedure
double get_objective_coefficient(const vector<vector<int>>& dist, const vector<int>& population, int i, int j);

//find dual vars of hess model for Lagrangian hot start
hess_params findHotDual(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1);

// build hess model and return x variables
hess_params build_hess(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1);
// add MCF constraints to model with hess variables x
void build_shir(GRBModel* model, hess_params& p, graph* g);
void build_mcf(GRBModel* model, hess_params& p, graph* g);
// add CUT constraints to model with hess variables x (lazy)
class HessCallback : public GRBCallback
{
protected:
  hess_params& p;
  double** x_val; // x values
  graph* g; // graph pointer
  int n; // g->nr_nodes
  const vector<int> population;
public:
  int numCallbacks; // number of callback calls
  double callbackTime; // cumulative time in callbacks
  int numLazyCuts;
  HessCallback(hess_params& p_, graph* g_, const vector<int>& population_) : p(p_), g(g_), population(population_), numCallbacks(0), callbackTime(0.), numLazyCuts(0)

  {
    n = g->nr_nodes;
    x_val = new double*[n];
    for (int i = 0; i < n; ++i)
      x_val[i] = new double[n];
  }
  virtual ~HessCallback()
  {
    for (int i = 0; i < n; ++i)
      delete[] x_val[i];
    delete[] x_val;
  }
protected:
  void populate_x()
  {
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (IS_X(i, j))
          x_val[i][j] = getSolution(X_V(i, j));
        else if (p.F1[i][j]) x_val[i][j] = 1.; else x_val[i][j] = 0.;
  }
};

// @return callback for delete only
HessCallback* build_cut(GRBModel* model, hess_params& p, graph* g, const vector<int>& population);
HessCallback* build_lcut(GRBModel* model, hess_params& p, graph* g, const vector<int>& population, int U);
//Lagrangian functions
// input:
//    g: graph pointer
//    multipliers : 3|V| array of lagrangian multipliers [A,L,U]
//    L : population lower limit
//    U : population upper limit
//    k : number of districts
//    population : array i-th element is population of i-th node
//    w : original objective coefficients w_ij to assign i to j
//    w_hat : adjusted objective coefficients (after combining like terms)
//    W : weight of a min-weight subgraph rooted at j is W_j, i.e., w_hat_jj + \sum_{j!=i} max(0,w_hat_ij)
//    S : solution of the inner problem for clusterheads
//    grad : pointer to the resulting gradient
//    f_val : resulting objective value
//    currentCenters : the best k centers (for the current multipliers), i.e., the k vertices j that have least W_j
void solveInnerProblem(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
  const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& currentCenters);

double solveLagrangian(graph* g, const vector<vector<double>>& w, const vector<int> &population, int L, int U, int k,
  vector<vector<double>>& LB1, bool ralg_hot_start, const char* ralg_hot_start_fname, const run_params& rp, bool exploit_contiguity);

void update_LB(const vector<double>& W, const vector<bool>& currentCenters, double f_val,
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB1);

void update_LB_contiguity(graph* g, const vector<double>& W, const vector<bool>& currentCenters, double f_val,
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB1);

vector<int> HessHeuristic(graph* g, const vector<vector<double> >& w, const vector<int>& population,
  int L, int U, int k, double &UB, int maxIterations, bool do_cuts = false);

void ContiguityHeuristic(vector<int> &heuristicSolution, graph* g, const vector<vector<double> > &w, 
  const vector<int> &population, int L, int U, int k, double &UB, string arg_model);

bool LocalSearch(graph* g, const vector<vector<double> >& w, const vector<int>& population,
  int L, int U, int k, vector<int>&heuristicSolution, double &UB);

#endif
