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
void build_mcf0(GRBModel* model, GRBVar** x, graph* g);
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
        for (int i = 0; i < n; ++i)
            x[i] = new double[n];
    }
    virtual ~HessCallback()
    {
        for (int i = 0; i < n; ++i)
            delete[] x[i];
        delete[] x;
    }
protected:
    void populate_x()
    {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                x[i][j] = getSolution(grb_x[i][j]);
    }
};

// @return callback for delete only
HessCallback* build_cut1(GRBModel* model, GRBVar** x, graph* g);
HessCallback* build_cut2(GRBModel* model, GRBVar** x, graph* g);

// add UL1 & UL2 instance and return x variables
GRBVar** build_UL_1(GRBModel* model, graph* g, const vector<int>& population, int k);
GRBVar** build_UL_2(GRBModel* model, graph* g, const vector<int>& population, int k);

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
void solveInnerProblem(graph* g, const double* multipliers, const vector<vector<bool>>& F_0, const vector<vector<bool>>& F_1,
    int L, int U, int k, const vector<vector<int>>& clusters, const vector<int>& population,
    const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, vector<bool>& S, double* grad, double& f_val);

//Preprocess functions
vector<vector<int>> preprocess(graph* g, vector<int>& new_population, int L, int U, const vector<int>& population);
int FindMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, int L);
vector<vector<int>> FindClustersFromStemVector(graph* g, vector<int>& stem);
void QuickTestForInfeasibility(graph* g, vector<int>& new_population, vector<bool>& deleted, int L, int U);
void strengthen_hess(GRBModel* model, GRBVar** x, graph* g, vector<vector<int>>& clusters);
#endif
