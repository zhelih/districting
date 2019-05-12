#ifndef _MODELS_H
#define _MODELS_H

#include <vector>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

typedef const vector<vector<bool>> cvv;

//hack
#define IS_X(i,j) (!F0[i][j] && !F1[i][j])
#define X(i,j) (F0[i][j]?GRBLinExpr(0.):(F1[i][j]?GRBLinExpr(1.):GRBLinExpr(x[i][j])))

//auxilliary procedure
double get_objective_coefficient(const vector<vector<int>>& dist, const vector<int>& population, int i, int j);

// build hess model and return x variables
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1);
// add SCF constraints to model with hess variables x
void build_scf(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);
// add MCF constraints to model with hess variables x
void build_mcf0(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);
void build_mcf1(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);
void build_mcf2(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);
// add CUT constraints to model with hess variables x (lazy)
class HessCallback : public GRBCallback
{
protected:
    GRBVar** x; // x variables
    double** x_val; // x values
    graph* g; // graph pointer
    int n; // g->nr_nodes
    cvv& F0;
    cvv& F1;
public:
    int numCallbacks; // number of callback calls
    double callbackTime; // cumulative time in callbacks
    int numLazyCuts;
    HessCallback(GRBVar** grb_x_, graph* g_, cvv& F0_, cvv& F1_) : x(grb_x_), g(g_), F0(F0_), F1(F1_), numCallbacks(0), callbackTime(0.), numLazyCuts(0)
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
              if(IS_X(i,j))
                x_val[i][j] = getSolution(x[i][j]);
              else if(F1[i][j]) x_val[i][j] = 1.; else x_val[i][j] = 0.;
    }
};

// @return callback for delete only
HessCallback* build_cut1(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);
HessCallback* build_cut2(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1);

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
void eugene_inner(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
    const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& S,
    const vector<vector<bool>>& F_0, const vector<vector<bool>>& F_1);

void lagrangianBasedSafeFixing(vector<vector<bool>>& F_0, vector<vector<bool>>& F_1,
    const vector<vector<int>>& clusters, vector<double>& W, const vector<bool>& S, const double f_val, const double UB, const vector<vector<double>> &w_hat);

//Preprocess functions
vector<vector<int>> preprocess(graph* g, vector<int>& new_population, int L, int U, const vector<int>& population);
int FindMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, int L);
vector<vector<int>> FindClustersFromStemVector(graph* g, vector<int>& stem);
void QuickTestForInfeasibility(graph* g, vector<int>& new_population, vector<bool>& deleted, int L, int U);
void strengthen_hess(GRBModel* model, GRBVar** x, graph* g, vector<vector<int>>& clusters);


#endif
