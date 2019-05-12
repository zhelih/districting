// source file for building hess model for further usage
#include <cstdio>
#include <vector>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"

using namespace std;

double get_objective_coefficient(const vector<vector<int>>& dist, const vector<int>& population, int i, int j)
{
  return (static_cast<double>(dist[i][j])/1000.) * (static_cast<double>(dist[i][j])/1000.) * static_cast<double>(population[i]);
}


// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1)
{
  // create GUROBI Hess model
  int n = g->nr_nodes;

  // create n^2 variables
  GRBVar** x = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
  {
    x[i] = new GRBVar[n];
    for(int j = 0; j < n; ++j)
      if(IS_X(i,j))
        x[i][j] = model->addVar(0., 1., 0., GRB_BINARY);
  }
  model->update();

  // Set objective: minimize sum d^2_ij*x_ij
  GRBLinExpr expr = 0;
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      expr += get_objective_coefficient(dist, population, i, j) * X(i,j);

  model->setObjective(expr, GRB_MINIMIZE);

  // add constraints (1b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for(int j = 0; j < n; ++j)
        constr += X(i,j);
    model->addConstr(constr == 1);
  }

  // add constraint (1c)
  expr = 0;
  for(int j = 0; j < n; ++j)
    expr += X(j,j);
  model->addConstr(expr == k);

  // add constraint (1d)
  for(int j = 0; j < n; ++j)
  {
    GRBLinExpr constr = 0;
    for(int i = 0; i < n; ++i)
      constr += population[i]*X(i,j);
    // add for j
    model->addConstr(constr - U*X(j,j) <= 0); // U
    model->addConstr(constr - L*X(j,j) >= 0); // L
  }

  // add contraints (1e)
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      if(i != j)
        model->addConstr(X(i,j) <= X(j,j));

  model->update();

  model->write("debug_hess.lp");

  return x;
}

