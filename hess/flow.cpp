// source file for single and multi commodity flow formulations
#include <unordered_map>
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

  // step 1 : hash edge (i,j) to i*n+j = h1
  // step 2 : hash every h1 to a number, resulting in exactly |E| variables

  std::unordered_map<int,int> hash_edges;
  int cur = 0;
  for(int i = 0; i < n; ++i)
    for(int j : g->nb(i))
      hash_edges.insert(std::make_pair(i*n+j, cur++));

  // debug
/*  for(int i = 0; i < n; ++i)
    for(int j : g->nb(i))
      printf("edge (%d,%d) hashed as %d\n", i, j, hash_edges[i*n+j]);
*/

  // get number of edges
  int nr_edges = hash_edges.size();

  // add additional flow variable (many) f[v][i,j]
  GRBVar**f = new GRBVar*[n]; // commodity type, v
  for(int v = 0; v < n; ++v)
    f[v] = model->addVars(nr_edges, GRB_CONTINUOUS); // the edge

  // update variables to binary when j == v
  for(int v = 0; v < n; ++v)
    for(int i = 0; i < n; ++i)
      for(int j : g->nb(i))
        if(v == j)
          f[v][hash_edges[n*i+j]].set(GRB_CharAttr_VType, GRB_BINARY);

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
        expr += f[i][hash_edges[n*j+nb_j]]; // in d_+ : edge (j -- nb_j)
        expr -= f[i][hash_edges[n*nb_j+j]]; // in d_- : edge (nb_j -- j)
      }
      model->addConstr(expr - x[i][j] == 0);
    }
  }

  // add constraint (16c)
  for(int i = 0; i < n; ++i)
  {
    GRBLinExpr expr = 0;
    for(int nb_i : g->nb(i))
      expr += f[i][hash_edges[n*i+nb_i]]; // in d_+ : edge (i -- nb_i)
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
        model->addConstr(f[v][hash_edges[n*i+j]] - f[j][hash_edges[n*i+j]] <= 0);
      }

  model->update();

  model->write("debug_mcf1.lp");
}

void build_mcf2(GRBModel* model, GRBVar** x, graph* g)
{
  throw "Unimplemented!";
}
