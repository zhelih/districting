#ifndef _MODELS_H
#define _MODELS_H

#include <vector>
#include <unordered_map>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

typedef const vector<vector<bool>> cvv;

struct hess_params
{
  GRBVar* x;
  vector<vector<bool>> F0;
  vector<vector<bool>> F1;
  unordered_map<int, int> h;
  int n;
};

//hack
#define IS_X(i,j) (!p.F0[i][j] && !p.F1[i][j])
#define X_V(i,j) (p.x[p.h[p.n*i+j]])
#define X(i,j) (p.F0[i][j]?GRBLinExpr(0.):(p.F1[i][j]?GRBLinExpr(1.):GRBLinExpr(X_V(i,j))))

//auxilliary procedure
double get_objective_coefficient(const vector<vector<int>>& dist, const vector<int>& population, int i, int j);

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
//    F_0 : vertices fixed x_jj = 0
//    F_1 : vertices fixed x_jj = 1
//    L : population lower limit
//    U : population upper limit
//    k : number of districts
//    cluters : ?
//    population : array i-th element is population of i-th node
//    w : ?
//    w_hat : ?
//    W : ?
//    S : solution of the inner problem for clusterheads
//    grad : pointer to the resulting gradient
//    f_val : resulting objective value
void eugene_inner(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
  const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& S,
  const vector<vector<bool>>& F_0, const vector<vector<bool>>& F_1);

void lagrangianBasedSafeFixing(vector<vector<bool>>& F_0, vector<vector<bool>>& F_1,
  const vector<vector<int>>& clusters, vector<double>& W, const vector<bool>& S, const double f_val, const double UB, const vector<vector<double>> &w_hat);

double solveLagrangian(graph* g, const vector<vector<double>>& w, const vector<int> &population, int L, int U, int k,
  vector<vector<double>>& LB0, vector<vector<double>>& LB1, vector<int> &lagrangianCenters, bool ralg_hot_start, const char* ralg_hot_start_fname, bool exploit_contiguity);

void solveInnerProblem(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
  const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& currentCenters);

void update_LB0_and_LB1(const vector<double>& W, const vector<bool>& currentCenters, double f_val,
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB0, vector< vector<double> > &LB1);

void update_LB0_and_LB1_contiguity(graph* g, const vector<double>& W, const vector<bool>& currentCenters, double f_val,
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB0, vector< vector<double> > &LB1);

vector<int> HessHeuristic(graph* g, const vector<vector<double> >& w, const vector<int>& population,
  int L, int U, int k, double &UB, int maxIterations, bool do_cuts = false);

void ContiguityHeuristic(vector<int> &heuristicSolution, graph* g, const vector<vector<double> > &w, 
  const vector<int> &population, int L, int U, int k, double &UB, string arg_model);

void LocalSearch(graph* g, const vector<vector<double> >& w, const vector<int>& population,
  int L, int U, int k, vector<int>&heuristicSolution, string arg_model, double &UB);// , cvv &F0);

//Preprocess functions
vector<vector<int>> preprocess(graph* g, vector<int>& new_population, int L, int U, const vector<int>& population);
int FindMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, int L);
vector<vector<int>> FindClustersFromStemVector(graph* g, vector<int>& stem);
void QuickTestForInfeasibility(graph* g, vector<int>& new_population, vector<bool>& deleted, int L, int U);
void strengthen_hess(GRBModel* model, hess_params& p, graph* g, vector<vector<int>>& clusters);

#endif
