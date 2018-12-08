// source file for single and multi commodity flow formulations
#include "gurobi_c++.h"
#include "graph.h"

void build_scf(GRBModel* model, GRBVar** x, graph* g)
{
  // create n^2 variables for arcs, presolve will eliminate unused
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
  for(int i = 0; i < n; ++i)
    for(int j : g->nb(i))
    {
      model->addConstr(f[i][j] - y[i][j] >= 0);
      model->addConstr(f[i][j] - n*y[i][j] <= 0);
    }

  // add constraint (Austin)
  for(int v = 0; v < n; ++v)
  {
    // for every edge here
    for(int i = 0; i < n; ++i)
      for(int j : g->nb(i))
        model->addConstr(x[i][v] + y[i][j] - x[j][v] <= 1);
  }
  model->update();

  model->write("debug_scf.lp");
}

void build_mcf1(GRBModel* model, GRBVar** x, graph* g)
{
  int n = g->nr_nodes;
  // add additional flow variable (many) f[v][i][j]
  GRBVar***f = new GRBVar**[n]; // commodity type, v
  for(int v = 0; v < n; ++v)
  {
    f[v] = new GRBVar*[n]; // from node, i
    for(int i = 0; i < n; ++i)
      f[v][i] = model->addVars(n, GRB_BINARY); // to node, j
  }

  model->update();

  // add constraint (16b)
  for(int i = 0; i < n; ++i)
  {
    for(int j = 0; j < n; ++j)
    {
      if(i == j)
        continue;
      GRBLinExpr expr = 0;
      for(int nb_j : g->nb(j))
      {
        expr += f[i][j][nb_j]; // in d_+ : edge (j -- nb_j)
        expr -= f[i][nb_j][j]; // in d_- : edge (nb_j -- j)
      }
      model->addConstr(expr - x[i][j] == 0);
    }
  }

  // add constraint (16c)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr expr = 0;
    for(int nb_i : g->nb(i))
      expr += f[i][i][nb_i]; // in d_+ : edge (i -- nb_i)
    model->addConstr(expr == 0);
  }

  // add constraint (16d)
  // FIXME: i < j only or nah?
  for(int i = 0; i < n; ++i)
    for(int j : g->nb(i))
      for(int v = 0; v < n; ++v)
      {
        if(v == i || v == j)
          continue;
        model->addConstr(f[v][i][j] - f[j][i][j] <= 0);
      }

  model->update();

  model->write("debug_mcf1.lp");
}

void build_mcf2(GRBModel* model, GRBVar** x, graph* g)
{
  throw "Unimplemented!";
}
