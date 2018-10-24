// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

int read_input_data(const char* dimacs_fname, const char* distance_fname, const char* population_fname, // INPUTS
                     graph* g, vector<vector<int> >& dist, vector<int> population) // OUTPUTS
{
    // read dimacs graph
    g = from_dimacs(dimacs_fname);
    if(!g) {
      fprintf(stderr, "Failed to read dimacs graph from %s\n", dimacs_fname);
      return 1;
    }

    // read distances (must be sorted)
    FILE* f = fopen(distance_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", distance_fname);
      return 1;
    }
    // file contains the first row as a header row, skip it
    // also skip the first element in each row (node id)
    char buf[3000]; //dummy
    fgets(buf, sizeof(buf), f); // skip first line
    dist.resize(g->nr_nodes);
    for(int i = 0; i < g->nr_nodes; ++i)
    {
      dist[i].resize(g->nr_nodes);
      int d; fscanf(f, "%d,", &d); // skip first element
      for(int j = 0; j < g->nr_nodes; ++j) {
        fscanf(f, "%d,", &d);
        dist[i][j] = d;
      }
    }
    fclose(f);

    // read population file
    f = fopen(population_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", population_fname);
      return 1;
    }
    // skip first line about total population
    fgets(buf, sizeof(buf), f);
    population.resize(g->nr_nodes);
    for(int i = 0; i < g->nr_nodes; ++i) {
      int node, pop;
      fscanf(f, "%d %d ", &node, &pop);
      population[node] = pop;
    }
    fclose(f);
    return 0;
}

void run_gurobi(graph* g, const vector<vector<int> >& dist, const vector<int>& population, int L, int U, int k)
{
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
        expr += ((double)dist[i][j]/1000.) * ((double)dist[i][j]/1000.) * population[i] * x[i][j]; //FIXME check overflow?

    model.setObjective(expr, GRB_MINIMIZE);

    // add constraints (1b)
    for(int i = 0; i < n; ++i)
    {
      GRBLinExpr constr = 0;
      for(int j = 0; j < n; ++j)
          constr += x[i][j];
      model.addConstr(constr == 1);
    }

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

    // add contraints (1e)
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(i != j)
          model.addConstr(x[i][j] <= x[j][j]);

    // Optimize model
    model.write("debug.lp");
    model.optimize();

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

    // read inputs
    graph* g = 0;
    vector<vector<int> > dist;
    vector<int> population;
    if(read_input_data(dimacs_fname, distance_fname, population_fname, g, dist, population))
      return 1; // failure

    // check connectivity
    if(!g->is_connected())
    {
      printf("Problem is infeasible (not connected!)\n");
      return 0;
    }

    run_gurobi(g, dist, population, L, U, k);

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
