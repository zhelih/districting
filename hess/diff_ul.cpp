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
    p_bar += static_cast<double>(p) / static_cast<double>(k);

  // create n^2 variables
  GRBVar** x = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    x[i] = model->addVars(n, GRB_BINARY);

  GRBVar L = *model->addVars(1, GRB_CONTINUOUS);
  GRBVar U = *model->addVars(1, GRB_CONTINUOUS);
  model->update();

  model->setObjective(U-L, GRB_MINIMIZE);

  // add constraints (36b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for(int j = 0; j <= i; ++j)
        constr += x[i][j];
    model->addConstr(constr == 1);
  }

  // add constraint (36c)
  GRBLinExpr expr = 0;
  for(int j = 0; j < n; ++j)
    expr += x[j][j];
  model->addConstr(expr == k);

  // add constraint (36d)
  for(int j = 0; j < n; ++j)
  {
    GRBLinExpr constr = 0;
    for(int i = j; i < n; ++i)
      constr += population[i]*x[i][j];
    // add for j
    model->addConstr(constr - U + p_bar * x[j][j] <= p_bar);
    model->addConstr(constr - L - p_bar * x[j][j] >= -p_bar);
  }

  // add contraints (36e)
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < i; ++j)
      model->addConstr(x[i][j] <= x[j][j]);

  // add constraint (36f)
  model->addConstr(L <= p_bar);
  // add constraint (36g)
  model->addConstr(U >= p_bar);

  model->update();
  model->write("debug_ul1.lp");

  return x;
}

// add UL2 instance
GRBVar** build_UL_2(GRBModel* model, graph* g, const vector<int>& population, int k)
{
  int n = g->nr_nodes;

  GRBVar** x = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    x[i] = model->addVars(k, GRB_BINARY);

  GRBVar* lu = model->addVars(2, GRB_CONTINUOUS); // L and U as variables
  model->update();

  model->setObjective(lu[1] - lu[0], GRB_MINIMIZE);

  // add constraints (37b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for(int j = 0; j < k; ++j)
        constr += x[i][j];
    model->addConstr(constr == 1);
  }

  // add constraint (37c)
  for(int j = 0; j < k; ++j)
  {
    GRBLinExpr constr = 0;
    for(int i = 0; i < n; ++i)
      constr += population[i]*x[i][j];
    // add for j
    model->addConstr(constr - lu[1] <= 0); // U
    model->addConstr(constr - lu[0] >= 0); // L
  }

  // explicitly fix for zero
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < k; ++j)
      if(i < j)
        model->addConstr(x[i][j] == 0);

  model->set(GRB_IntParam_Symmetry, 2);

  model->update();
  model->write("debug_ul2.lp");

  return x;
}
