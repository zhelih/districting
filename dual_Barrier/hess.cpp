// source file for building hess model for further usage
#include <cstdio>
#include <vector>
#include "graph.h"
#include "gurobi_c++.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::sort


using namespace std;

// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k)
{
    // create GUROBI Hess model
    int n = g->nr_nodes;

    // create n^2 variables
    GRBVar** x = new GRBVar*[n];
    for (int i = 0; i < n; ++i)
        x[i] = model->addVars(n, GRB_CONTINUOUS);
    model->update();

    // Set objective: minimize sum d^2_ij*x_ij
    GRBLinExpr expr = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            expr += ((double)dist[i][j] / 1000.) * ((double)dist[i][j] / 1000.) * population[i] * x[i][j]; //FIXME check overflow?

    model->setObjective(expr, GRB_MINIMIZE);

    //x[452][452].set(GRB_DoubleAttr_LB, 1.0);
    //x[429][429].set(GRB_DoubleAttr_LB, 1.0);

    // add constraints (1b)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr constr = 0;
        for (int j = 0; j < n; ++j)
            constr += x[i][j];
        model->addConstr(constr == 1);
    }

 

    // add constraint (1d-LB)
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * x[i][j];
        // add for j
        //model->addConstr(constr - U * x[j][j] <= 0); // U
        model->addConstr(constr - L * x[j][j] >= 0); // L
    }

    // add constraint (1d-UB)
    for (int j = 0; j < n; ++j)
    {
        GRBLinExpr constr = 0;
        for (int i = 0; i < n; ++i)
            constr += population[i] * x[i][j];
        // add for j
        model->addConstr(constr - U * x[j][j] <= 0); // U
        //model->addConstr(constr - L * x[j][j] >= 0); // L
    }

    // add constraint (1c)
    expr = 0;
    for (int j = 0; j < n; ++j)
        expr += x[j][j];
    model->addConstr(expr == k);

    // add contraints (1e)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                model->addConstr(x[i][j] <= x[j][j]);

    model->update();

    model->write("debug_hess.lp");

    return x;
}

string statusNumtoString(int num)
{
    if (num == 1) return "Model is loaded, but no solution information is available.";
    else if (num == 2) return "Model was solved to optimality (subject to tolerances), and an optimal solution is available.";
    else if (num == 3) return "Model was proven to be infeasible.";
    else if (num == 4) return "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.";
    else if (num == 5) return "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.";
    else if (num == 6) return "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.";
    else if (num == 7) return "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.";
    else if (num == 8) return "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.";
    else if (num == 9) return "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.";
    else if (num == 10) return "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.";
    else if (num == 11) return "Optimization was terminated by the user.";
    else if (num == 12) return "Optimization was terminated due to unrecoverable numerical difficulties.";
    else if (num == 13) return "Unable to satisfy optimality tolerances; a sub-optimal solution is available.";
    else if (num == 14) return "An asynchronous optimization call was made, but the associated optimization run is not yet complete.";
    else if (num >= 15 || num <= 0) return "No specific error could be recognized.";
}

