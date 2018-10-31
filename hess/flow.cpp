// source file for single and multi commodity flow formulations
#include "gurobi_c++.h"
#include "graph.h"

void build_scf(GRBModel* model, GRBVar** x, graph* g)
{
  // create n^2 variables for arcs, presolve will eliminate unused
  // TODO sophisticated indexing?
  int n = g->nr_nodes;

  GRBVar** y = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    y[i] = model->addVars(n, GRB_BINARY);
  GRBVar** f = new GRBVar*[n];
  for(int i = 0; i < n; ++i)
    f[i] = model->addVars(n, GRB_CONTINUOUS);

  model->update();

  // add constraint (16b)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr expr = 0;
    for(int j : g->nb(i))
    {
      expr += y[j][i];
    }
    model->addConstr(expr + x[i][i] == 1);
  }

  // add constraint (16c)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr expr = 0;
    for(int j : g->nb(i))
    {
      expr += f[i][j]; // in d_+
      expr -= f[j][i]; // in d_-
    }

    for(int j = 0; j < n; ++j)
      expr -= x[j][i];

    model->addConstr(expr == -1);
  }

  // add constraint (16d)
  // TODO without is_edge
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      if(g->is_edge(i,j))
      {
        model->addConstr(f[i][j] - y[i][j] >= 0);
        model->addConstr(f[i][j] - n*y[i][j] <= 0);
      }
  model->update();
}

void build_mcf(GRBModel* model, GRBVar** x, graph* g)
{
  throw "Unimplemented!";
}
