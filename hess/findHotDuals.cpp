// source file for finding hot dual
#include <cstdio>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include "io.h"

using namespace std;

hess_params findHotDual(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1)
{

    // create GUROBI Hess model
    int n = g->nr_nodes;
    hess_params p;
    p.n = n;
    p.F0 = F0; p.F1 = F1; // copy

    // hash variables
    int cur = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (!F0[i][j] && !F1[i][j])
                p.h[n*i + j] = cur++; //FIXME implicit reuse of the map (i,j) -> n*i+j

    printf("Build hess : created %lu variables\n", p.h.size());
    int nr_var = static_cast<int>(p.h.size());

    // create variables
    p.x = model->addVars(nr_var, GRB_CONTINUOUS);
    model->update();

    // Set objective: minimize sum d^2_ij*x_ij
    GRBLinExpr expr = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            expr += w[i][j] * X(i, j);

    model->setObjective(expr, GRB_MINIMIZE);

    // add constraints (b)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr constr = 0;
        for (int j = 0; j < n; ++j)
            constr += X(i, j);
        model->addConstr(constr == 1);
    }



    // add constraint (d) for population lower bound
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * X(i, j);
        model->addConstr(constr - L * X(j, j) >= 0);
    }

    // add constraint (d) for population upper bound
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * X(i, j);
        model->addConstr(constr - U * X(j, j) <= 0);
    }

    // add constraint (c)
    expr = 0;
    for (int j = 0; j < n; ++j)
        expr += X(j, j);
    model->addConstr(expr == k);

    // add contraints (e)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j && !F0[i][j])
                model->addConstr(X(i, j) <= X(j, j));

    model->update();

    return p;
}

