// source file for building hess model for further usage
#include <cstdio>
#include <vector>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k)
{
  // create GUROBI Hess model
  int n = g->nr_nodes;
  if(n <= 0)
  {
    fprintf(stderr, "build_hess: empty graph\n");
    throw "build_hess failed";
  }

  if(dist.size() != n || population.size() != n)
  {
    fprintf(stderr, "build_hess: dist/population size != n, expected %d\n", n);
    throw "build_hess failed";
  }

  // create n^2 variables
  GRBVar** x = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    x[i] = model->addVars(n, GRB_BINARY);
  model->update();

  // Set objective: minimize sum d^2_ij*x_ij
  GRBLinExpr expr = 0;
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      expr += ((double)dist[i][j]/1000.) * ((double)dist[i][j]/1000.) * population[i] * x[i][j]; //FIXME check overflow?

  model->setObjective(expr, GRB_MINIMIZE);

  // add constraints (1b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for(int j = 0; j < n; ++j)
        constr += x[i][j];
    model->addConstr(constr == 1);
  }

  // add constraint (1c)
  expr = 0;
  for(int j = 0; j < n; ++j)
    expr += x[j][j];
  model->addConstr(expr == k);

  // add constraint (1d)
  for(int j = 0; j < n; ++j)
  {
    GRBLinExpr constr = 0;
    for(int i = 0; i < n; ++i)
      constr += population[i]*x[i][j];
    // add for j
    model->addConstr(constr - U*x[j][j] <= 0); // U
    model->addConstr(constr - L*x[j][j] >= 0); // L
  }

  // add contraints (1e)
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      if(i != j)
        model->addConstr(x[i][j] <= x[j][j]);

  model->write("debug_hess.lp");

  return x;
}

