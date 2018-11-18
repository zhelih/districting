// source file for IP instances for best U-L difference
#include <cstdio>
#include <vector>
#include "graph.h"
#include "gurobi_c++.h"

using namespace std;

// add UL1 instance
// returns "x" variables
GRBVar** build_UL_1(GRBModel* model, graph* g, const vector<int>& population, int k)
{
  int n = g->nr_nodes;

  // compute p_bar
  double p_bar = 0.;
  for(int p: population)
    p_bar += p / k;

  // create n^2 variables
  GRBVar** x = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    x[i] = model->addVars(n, GRB_BINARY);

  GRBVar* lu = model->addVars(2, GRB_CONTINUOUS);
  model->update();

  model->setObjective(lu[1] - lu[0], GRB_MINIMIZE);

  // add constraints (37b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for(int j = 0; j <= i; ++j)
        constr += x[i][j];
    model->addConstr(constr == 1);
  }

  // add constraint (37c)
  GRBLinExpr expr = 0;
  for(int j = 0; j < n; ++j)
    expr += x[j][j];
  model->addConstr(expr == k);

  // add constraint (37d)
  for(int j = 0; j < n; ++j)
  {
    GRBLinExpr constr = 0;
    for(int i = j; i < n; ++i)
      constr += population[i]*x[i][j];
    // add for j
    model->addConstr(constr - lu[1] + p_bar * x[j][j] <= p_bar); // U
    model->addConstr(constr - lu[0] - p_bar * x[j][j] >= -p_bar); // L
  }

  // add contraints (37e)
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < i; ++j)
      model->addConstr(x[i][j] <= x[j][j]);

  // add constraint (37f)
  model->addConstr(lu[0] <= p_bar);
  // add constraint (37g)
  model->addConstr(lu[1] >= p_bar);

  model->update();
  model->write("debug_lu1.lp");

  return x;
}
