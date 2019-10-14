#include "gurobi_c++.h"
#include "graph.h"
#include "io.h"
#include "models.h"
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>
#include <string>
#include "common.h"

using namespace std;
extern const char* gitversion;

int main(int argc, char *argv[])
{
  printf("Districting, build %s\n", gitversion);
  if (argc < 2) {
    printf("Usage: %s <config> [state [ralg_hot_start]]\n", argv[0]);
    return 0;
  }

  // parse config
  run_params rp;
  rp = read_config(argv[1], (argc>2 ? argv[2] : ""), (argc>3 ? argv[3] : ""));
  int L = rp.L; int U = rp.U;
  const char* ralg_hot_start_fname = (rp.ralg_hot_start.empty() ? nullptr : rp.ralg_hot_start.c_str());

  if(!ralg_hot_start_fname)
  {
    printf("I have no idea where to dump ralg_hot_start, refusing to run!\n");
    return 1; // fail
  }

  // read inputs
  graph* g = nullptr;
  vector<vector<int> > dist;
  vector<int> population;
  if (read_input_data(rp.dimacs_file.c_str(), rp.distance_file.c_str(), rp.population_file.c_str(), g, dist, population))
    return 1; // fail

  int k = (rp.k == 0) ? g->get_k() : rp.k;

  if (L == 0 || U == 0)
    calculate_UL(population, k, &L, &U);

  printf("Model input: L = %d, U = %d, k = %d.\n", L, U, k);
  g->connect(dist);

  int nr_nodes = g->nr_nodes;
  // set objective function coefficients
  vector<vector<double>> w(nr_nodes, vector<double>(nr_nodes)); // this is the weight matrix in the objective function
  for (int i = 0; i < nr_nodes; i++)
    for (int j = 0; j < nr_nodes; j++)
      w[i][j] = get_objective_coefficient(dist, population, i, j);


  // determine which variables can be fixed
  vector<vector<bool>> F0(nr_nodes, vector<bool>(nr_nodes, false)); // define matrix F_0
  vector<vector<bool>> F1(nr_nodes, vector<bool>(nr_nodes, false)); // define matrix F_1

  auto start = chrono::steady_clock::now();

  try
  {
    // initialize environment and create an empty model
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // get incumbent solution using centers from lagrangian
    hess_params p;
    p = build_hess_special(&model, g, w, population, L, U, k); // constraints are well-organized

    // relax the model
    //model.set(GRB_IntParam_NodeCount, 1);
    for(int i = 0; i < nr_nodes; ++i)
      for(int j = 0; j < nr_nodes; ++j)
      {
        X_V(i, j).set(GRB_CharAttr_VType, 'C'); // set variables continuous
        X_V(i, j).set(GRB_DoubleAttr_UB, 1.); // upper bound
      }

    model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour
    //model.set(GRB_IntParam_Threads, 10); // limit to 10 threads
    model.set(GRB_DoubleParam_NodefileStart, 10); // 10 GB
    model.set(GRB_IntParam_Method, 3);  // use concurrent method to solve root LP

    model.optimize();

    double opt = model.get(GRB_DoubleAttr_ObjVal);
    GRBConstr* c = model.getConstrs();
    // extra memory usage, but should not matter
    std::vector<double> x_val(3*nr_nodes);
    for(int i = 0; i < 3*nr_nodes; ++i)
    {
      double coef = 1;
      if(i >= nr_nodes) coef = L;
      if(i >= 2*nr_nodes) coef = U;
      x_val[i] = coef * c[i].get(GRB_DoubleAttr_Pi);
    }
    dump_ralg_hot_start_fname(ralg_hot_start_fname, x_val.data(), 3*nr_nodes, opt);

    chrono::duration<double> duration = chrono::steady_clock::now() - start;
    printf("Total time elapsed: %lf seconds\n", duration.count());
  }
  catch (GRBException e) {
    printf("Error code = %d\n", e.getErrorCode());
    printf("%s\n", e.getMessage().c_str());
  }
  catch (const char* msg) {
    printf("Exception with message : %s\n", msg);
  }
  catch (...) {
    printf("Exception during optimization\n");
  }

  if(g) delete g;
  return 0;
}
