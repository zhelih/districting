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

void find_LP_dual_solution(GRBModel* model, graph* g, run_params rp, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k)
{

    // create GUROBI Hess model
    int n = g->nr_nodes;
    hess_params p;
    p.n = n;
    //p.F0 = F0; p.F1 = F1; // copy

    // hash variables
    int cur = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            //cerr << i << ", " << j << endl;
            //if (!F0[i][j] && !F1[i][j])
            p.h[n*i + j] = cur++; //FIXME implicit reuse of the map (i,j) -> n*i+j
        }            
    }
        

    printf("Build hess : created %lu variables\n", p.h.size());
    int nr_var = static_cast<int>(p.h.size());

    // create variables
    p.x = model->addVars(nr_var, GRB_CONTINUOUS);
    model->update();

    //cerr << "Hi there!" << endl;

    // Set objective: minimize sum d^2_ij*x_ij
    GRBLinExpr expr = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            //cerr << i << ", " << j << endl;
            expr += w[i][j] * X_V(i, j);
        }        
    }
        

    model->setObjective(expr, GRB_MINIMIZE);

    // add constraints (b)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr constr = 0;
        for (int j = 0; j < n; ++j)
            constr += X_V(i, j);
        model->addConstr(constr == 1);
    }



    // add constraint (d) for population lower bound
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * X_V(i, j);
        model->addConstr(constr - L * X_V(j, j) >= 0);
    }

    // add constraint (d) for population upper bound
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * X_V(i, j);
        model->addConstr(constr - U * X_V(j, j) <= 0);
    }

    // add constraint (c)
    expr = 0;
    for (int j = 0; j < n; ++j)
        expr += X_V(j, j);
    model->addConstr(expr == k);

    // add contraints (e)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            //if (i != j && !F0[i][j])
                model->addConstr(X_V(i, j) <= X_V(j, j));

    model->update();

    model->set(GRB_IntParam_Method, 2);
    model->set(GRB_IntParam_OutputFlag, 0);
    model->set(GRB_IntParam_Crossover, 0);

    model->optimize();

    if (model->get(GRB_IntAttr_Status) == 3) // infeasible
        cerr << "Infeasible!" << endl;

    //define constraints
    GRBConstr *c = 0;
    c = model->getConstrs();

    int i;

    //create a file
    FILE* f;
    string dualHot_fn = string(rp.state) + "_" + "hotDual" + ".hot";
    f = fopen(dualHot_fn.c_str(), "w");

    //print dual vars in file

    for (i = 0; i < n; ++i)
        fprintf(f, "%.6lf\n", c[i].get(GRB_DoubleAttr_Pi));

    for (; i < 2 * n; ++i)
        fprintf(f, "%.6lf\n", L * c[i].get(GRB_DoubleAttr_Pi));

    for (; i < 3 * n; ++i)
        fprintf(f, "%.6lf\n", U * c[i].get(GRB_DoubleAttr_Pi));

}
