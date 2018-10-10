// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

// obsolete for assymetric formulation where only i < j
// only for j < i
inline int ind(const int i, const int j, const int n)
{
  return j + (2*n - i - 1)*i/2;
}

//for debug
void print_GRBLinExpr(GRBLinExpr& expr)
{
  for(unsigned int i = 0; i < expr.size(); ++i)
  {
    if(expr.getCoeff(i) != 0)
      printf("+ %.1lf %s", expr.getCoeff(i), expr.getVar(i).get(GRB_StringAttr_VarName).c_str());
  }
  printf("\n");
//  printf(" <= %.1lf", expr.getConstant());
}

int main(int argc, char *argv[])
{
  ios::sync_with_stdio(1);
  if(argc < 3) {
    printf("Usage: %s <dimacs> <distance> <population> <L> <U> <k>\n", argv[0]);
    return 0;
  }
  try {
    // parse command line arguments
    char* dimacs_fname = argv[1];
    char* distance_fname = argv[2];
    char* population_fname = argv[3];
    int L = std::stoi(argv[4]);
    int U = std::stoi(argv[5]);
    int k = std::stoi(argv[6]);
    printf("Model input: L = %d, U = %d, k = %d\n", L, U, k);

    // read dimacs graph
    graph* g = from_dimacs(dimacs_fname);
    if(!g) {
      fprintf(stderr, "Failed to read dimacs graph from %s\n", dimacs_fname);
      return 1;
    }

    // read distances (must be sorted)
    vector<vector<int> > dist(g->nr_nodes, vector<int>(g->nr_nodes, 0));
    FILE* f = fopen(distance_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", distance_fname);
      return 1;
    }
    // file contains the first row as a header row, skip it
    // also skip the first element in each row (node id)
    char buf[3000]; //dummy
    fgets(buf, sizeof(buf), f); // skip first line
    for(int i = 0; i < g->nr_nodes; ++i)
    {
      int d; fscanf(f, "%d,", &d); // skip first element
      for(int j = 0; j < g->nr_nodes; ++j) {
        fscanf(f, "%d,", &d);
        dist[i][j] = d;
      }
    }
    fclose(f);

/* TODO remove
    for(int i = 0; i < g->nr_nodes; ++i) {
      for(int j = 0; j < g->nr_nodes; ++j)
        fprintf(stderr, "%d ", dist[i][j]);
      fprintf(stderr, "\n");
    }
*/

    // read population file
    vector<int> population(g->nr_nodes, 0);
    f = fopen(population_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", population_fname);
      return 1;
    }
    // skip first line about total population
    fgets(buf, sizeof(buf), f);
    for(int i = 0; i < g->nr_nodes; ++i) {
      int node, pop;
      fscanf(f, "%d %d ", &node, &pop);
      population[node] = pop;
    }
    fclose(f);

    // check connectivity
    if(!g->is_connected())
    {
      printf("Problem is infeasible (not connected!)\n");
      return 0;
    }

    // create GUROBI Hess model
    int n = g->nr_nodes;

    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
//    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    // create n^2 variables
    GRBVar** x = new GRBVar*[n];
    for(int i = 0; i < n; ++i)
      x[i] = model.addVars(n, GRB_BINARY);
   model.update();

    // Set objective: minimize sum d^2_ij*x_ij
    GRBLinExpr expr = 0;
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        expr += (dist[i][j]/100) * (dist[i][j]/100) * population[i] * x[i][j]; //FIXME check overflow?

    model.setObjective(expr, GRB_MINIMIZE);

    // add constraints (1b)
    for(int i = 0; i < n; ++i)
    {
      GRBLinExpr constr = 0;
      for(int j = 0; j < n; ++j)
          constr += x[i][j];
      model.addConstr(constr == 1);
    }

    // add contraints (1e)
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(i != j)
          model.addConstr(x[i][j] <= x[j][j]);

    // add constraint (1c)
    expr = 0;
    for(int j = 0; j < n; ++j)
      expr += x[j][j];
    model.addConstr(expr == k);

    // add constraint (1d)
    for(int j = 0; j < n; ++j)
    {
      GRBLinExpr constr = 0;
      for(int i = 0; i < n; ++i)
        constr += population[i]*x[i][j];
      // add for j
      model.addConstr(constr - U*x[j][j] <= 0); // U
      model.addConstr(constr - L*x[j][j] >= 0); // L
    }

    // Optimize model
    model.write("debug.lp");
    model.optimize();
/*
    // translate solution
    for(int i = 0; i < n; ++i)
      if(x[i][i][0].get(GRB_DoubleAttr_X) > 0.5)
      {
        printf("Cluster head %d:", i);
        for(int j = 0; j < n; ++j)
        {
          if(j != i && x[j][i] && x[j][i][0].get(GRB_DoubleAttr_X) > 0.5)
          {
            printf(" %u", j);
          }
        }
        printf("\n");
      }
//      for(int j = i; j < n; ++j)
//        printf("x[%d][%d] = %lf\n", i, j, x[ind(i,j,n)].get(GRB_DoubleAttr_X));

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    cout << "Number of callbacks = " << cb.numCallbacks << endl;
    cout << "Total time in callbacks = " << cb.callbackTime << " secs " << endl;

    // write to file for an output
    vector<vector<int> > clusters(n);;
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(x[i][j] && x[i][j][0].get(GRB_DoubleAttr_X) > 0.5)
          clusters[j].push_back(i);*/

/*    printf("Districts:\n");
    for(i = 0; i < nr_districts; ++i)
    {
      printf("District %d: (y[i] = %d)", i, ((y[i].get(GRB_DoubleAttr_X)>0.5)?1:0));
      for(auto it = clusters[i].begin(); it != clusters[i].end(); ++it)
        printf(" %d", *it);
      printf("\n");
    }
*/
/*    FILE* f = fopen("districting.out", "w");
    int nr = 0;
    for(int i = 0; i < n; ++i)
    {
      for(int node : clusters[i])
        fprintf(f, "%d %d\n", node, nr);
      nr += (clusters[i].size() > 0);
    }

    fclose(f);
*/

    delete g;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
