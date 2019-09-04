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

const double VarFixingEpsilon = 0.00001;

using namespace std;
extern const char* gitversion;

#define DO_BATCH_OUTPUT

int main(int argc, char *argv[])
{
  printf("Districting, build %s\n", gitversion);
  if (argc < 2) {
    printf("Usage: %s <config> [state [ralg_hot_start]]\n\
  Available models:\n\
  \thess\t\tHess model\n\
  \tshir\t\tHess model with SHIR\n\
  \tmcf\t\tHess model with MCF\n\
  \tcut\t\tHess model with CUT\n\
  \tlcut\t\tHess model with LCUT\n", argv[0]);
    return 0;
  }

  // parse config
  run_params rp;
  rp = read_config(argv[1], (argc>2 ? argv[2] : ""), (argc>3 ? argv[3] : ""));
  int L = rp.L; int U = rp.U;
  bool ralg_hot_start = !rp.ralg_hot_start.empty();
  const char* ralg_hot_start_fname = (rp.ralg_hot_start.empty() ? nullptr : rp.ralg_hot_start.c_str());

  // read inputs
  graph* g = 0;
  vector<vector<int> > dist;
  vector<int> population;
  if (read_input_data(rp.dimacs_file.c_str(), rp.distance_file.c_str(), rp.population_file.c_str(), g, dist, population))
    return 1; // failure

  int k = (rp.k == 0) ? g->get_k() : rp.k;

  if (L == 0 || U == 0)
    calculate_UL(population, k, &L, &U);

  printf("Model input: L = %d, U = %d, k = %d.\n", L, U, k);

  // dump run args to output
  ffprintf(rp.output, "%s, %s, %d, %d, %d, %d, ", rp.state, rp.model.c_str(), g->nr_nodes, k, L, U);

  g->connect(dist);

  // check connectivity
  if (!g->is_connected())
  {
    printf("Problem is infeasible (not connected!)\n");
    ffprintf(rp.output, "disconnected\n");
    fclose(rp.output);
    return 1;
  }

  if (g->nr_nodes <= 0)
  {
    printf("empty graph\n");
    ffprintf(rp.output, "empty graph\n");
    fclose(rp.output);
    return 1;
  }

  if (dist.size() != g->nr_nodes || population.size() != g->nr_nodes)
  {
    printf("dist/population size != n, expected %d\n", g->nr_nodes);
    ffprintf(rp.output, "bad input data\n");
    fclose(rp.output);
    return 1;
  }

  long total_pop = 0L;
  for (int p : population)
    total_pop += p;
  printf("Model input: total population = %ld\n", total_pop);

  string arg_model = rp.model;

  bool exploit_contiguity = false;
  if (arg_model != "hess")
    exploit_contiguity = true;

  // set objective function coefficients
  vector<vector<double>> w(g->nr_nodes, vector<double>(g->nr_nodes)); // this is the weight matrix in the objective function
  for (int i = 0; i < g->nr_nodes; i++)
    for (int j = 0; j < g->nr_nodes; j++)
      w[i][j] = get_objective_coefficient(dist, population, i, j);

  auto start = chrono::steady_clock::now();

  // apply Lagrangian 
  vector< vector<double> > LB1(g->nr_nodes, vector<double>(g->nr_nodes, -MYINFINITY)); // LB1[i][j] is a lower bound on problem objective if we fix x[i][j] = 1
  auto lagrange_start = chrono::steady_clock::now();
  double LB = solveLagrangian(g, w, population, L, U, k, LB1, ralg_hot_start, ralg_hot_start_fname, rp, exploit_contiguity);// lower bound on problem objective, coming from lagrangian
  chrono::duration<double> lagrange_duration = chrono::steady_clock::now() - lagrange_start;
  ffprintf(rp.output, "%.2lf, %.2lf, ", LB, lagrange_duration.count());

  auto dump_maybe_inf = [&rp](double val) { if (fabs(val-MYINFINITY) <= 1.) ffprintf(rp.output, "infinity, "); else ffprintf(rp.output, "%.2lf, ", val); };

  // run a heuristic
  double UB = MYINFINITY;
  int maxIterations = 10;   // 10 iterations is often sufficient
  auto heuristic_start = chrono::steady_clock::now();
  vector<int> heuristicSolution = HessHeuristic(g, w, population, L, U, k, UB, maxIterations, false);
  chrono::duration<double> heuristic_duration = chrono::steady_clock::now() - heuristic_start;
  dump_maybe_inf(UB);
  ffprintf(rp.output, "%.2lf, ", heuristic_duration.count());
  printf("Best solution after %d of HessHeuristic is %.2lf\n", maxIterations, UB);

  // run local search
  auto LS_start = chrono::steady_clock::now();
  bool ls_ok;
  ls_ok = LocalSearch(g, w, population, L, U, k, heuristicSolution, UB);
  chrono::duration<double> LS_duration = chrono::steady_clock::now() - LS_start;
  dump_maybe_inf(UB);
  ffprintf(rp.output, "%.2lf, ", LS_duration.count());
  printf("Best solution after local search is %.2lf\n", UB);

  if (arg_model != "hess" && ls_ok)  // solve contiguity-constrained problem, restricted to centers from heuristicSolution
  {
    UB = MYINFINITY;
    auto contiguity_start = chrono::steady_clock::now();
    ContiguityHeuristic(heuristicSolution, g, w, population, L, U, k, UB, "shir"); // arg_model);
    chrono::duration<double> contiguity_duration = chrono::steady_clock::now() - contiguity_start;
    dump_maybe_inf(UB);
    ffprintf(rp.output, "%.2lf, ", contiguity_duration.count());
  } else ffprintf(rp.output, "n/a, n/a, ");

  // determine which variables can be fixed
  vector<vector<bool>> F0(g->nr_nodes, vector<bool>(g->nr_nodes, false)); // define matrix F_0
  vector<vector<bool>> F1(g->nr_nodes, vector<bool>(g->nr_nodes, false)); // define matrix F_1
  for (int i = 0; i < g->nr_nodes; ++i)
    for (int j = 0; j < g->nr_nodes; ++j)
      if (LB1[i][j] > UB + VarFixingEpsilon) F0[i][j] = true;

  // report the number of fixings
  int numFixedZero = 0;
  int numFixedOne = 0;
  int numUnfixed = 0;
  int numCentersLeft = 0;
  for (int i = 0; i < g->nr_nodes; ++i)
  {
    if (!F0[i][i]) numCentersLeft++;
    for (int j = 0; j < g->nr_nodes; ++j)
    {
      if (F0[i][j]) numFixedZero++;
      else if (F1[i][j]) numFixedOne++;
      else numUnfixed++;
    }
  }
  printf("\n");
  printf("Number of variables fixed to zero = %d\n", numFixedZero);
  printf("Number of variables fixed to one  = %d\n", numFixedOne);
  printf("Number of variables not fixed     = %d\n", numUnfixed);
  printf("Number of centers left            = %d\n", numCentersLeft);
  double perc_var_fixed = (double)(g->nr_nodes*g->nr_nodes - numUnfixed) / (g->nr_nodes*g->nr_nodes);
  printf("Percentage of vars fixed = %.2lf\n", perc_var_fixed);
  ffprintf(rp.output, "%.2lf, ", perc_var_fixed);

  try
  {
    // initialize environment and create an empty model
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // get incumbent solution using centers from lagrangian
    hess_params p;
    p = build_hess(&model, g, w, population, L, U, k, F0, F1);

    // push GUROBI to branch over clusterheads
    for (int i = 0; i < g->nr_nodes; ++i)
      if (IS_X(i, i))
        X_V(i, i).set(GRB_IntAttr_BranchPriority, 1);

    HessCallback* cb = 0;

    if (arg_model == "shir")
      build_shir(&model, p, g);
    else if (arg_model == "mcf")
      build_mcf(&model, p, g);
    else if (arg_model == "cut")
      cb = build_cut(&model, p, g, population);
    else if (arg_model == "lcut")
      cb = build_lcut(&model, p, g, population, U);
    else if (arg_model != "hess") {
      printf("ERROR: Unknown model : %s\n", arg_model.c_str());
      exit(1);
    }

    //TODO change user-interactive?
    model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour
    //model.set(GRB_IntParam_Threads, 10); // limit to 10 threads
    model.set(GRB_DoubleParam_NodefileStart, 10); // 10 GB
    model.set(GRB_IntParam_Method, 3);  // use concurrent method to solve root LP
    model.set(GRB_DoubleParam_MIPGap, 0);  // force gurobi to prove optimality

    // provide IP warm start 
    if(ls_ok)
      for (int i = 0; i < g->nr_nodes; ++i)
      {
        for (int j = 0; j < g->nr_nodes; ++j)
          if (IS_X(i, j))
            X_V(i, j).set(GRB_DoubleAttr_Start, 0);

        if (IS_X(i, heuristicSolution[i]))
          X_V(i, heuristicSolution[i]).set(GRB_DoubleAttr_Start, 1);
      }

    //optimize the model
    auto IP_start = chrono::steady_clock::now();
    model.optimize();
    chrono::duration<double> IP_duration = chrono::steady_clock::now() - IP_start;
    ffprintf(rp.output, "%.2lf, ", IP_duration.count());
    printf("IP duration time: %lf seconds\n", IP_duration.count());
    chrono::duration<double> duration = chrono::steady_clock::now() - start;
    printf("Total time elapsed: %lf seconds\n", duration.count());
    ffprintf(rp.output, "%.2lf, ", duration.count());
    if (cb)
    {
      printf("Number of callbacks: %d\n", cb->numCallbacks);
      printf("Time in callbacks: %lf seconds\n", cb->callbackTime);
      printf("Number of lazy constraints generated: %d\n", cb->numLazyCuts);
      ffprintf(rp.output, "%d, %.2lf, %d, ", cb->numCallbacks, cb->callbackTime, cb->numLazyCuts);
      delete cb;
    } else ffprintf(rp.output, "n/a, n/a, n/a, ");

    // output overtly
    int max_pv = population[0];
    for (int pv : population)
      max_pv = max(max_pv, pv);


    ffprintf(rp.output, "%.2lf, ", static_cast<double>(max_pv) / static_cast<double>(U));

    if (model.get(GRB_IntAttr_Status) == 3) // infeasible
      ffprintf(rp.output, "infeasible, , , ");
    else {

      double objbound = model.get(GRB_DoubleAttr_ObjBound);

      // no incumbent solution was found, these values do no make sense
      if (model.get(GRB_IntAttr_SolCount) == 0)
        ffprintf(rp.output, "?, ?, %.2lf, ", objbound);
      else
      {
        double objval = model.get(GRB_DoubleAttr_ObjVal);
        double mipgap = model.get(GRB_DoubleAttr_MIPGap)*100.;
        ffprintf(rp.output, "%.2lf, %.2lf, ", objval, mipgap);
      }

      ffprintf(rp.output, "%.2lf, ", objbound);
    }

    long nodecount = static_cast<long>(model.get(GRB_DoubleAttr_NodeCount));
    ffprintf(rp.output, "%ld, ", nodecount);

    if(model.get(GRB_IntAttr_SolCount) > 0)
      ffprintf(rp.output, "Y");
    else
      ffprintf(rp.output, "N");
    // printf("qwerky567: sol: L = %.4lf, U = %.4lf\n", x[g->nr_nodes][0].get(GRB_DoubleAttr_X), x[g->nr_nodes][1].get(GRB_DoubleAttr_X));
    //else
    // printf("qwerly567: sol: no incumbent solution found!\n");

    if (model.get(GRB_IntAttr_Status) != 3) {
      vector<int> sol;
      translate_solution(p, sol, g->nr_nodes);
      string soln_fn = string(rp.state) + "_" + arg_model + ".sol";
      int len = soln_fn.length() + 1;
      char* char_array = new char[len];
      strcpy(char_array, soln_fn.c_str());
      printf_solution(sol, char_array);
    }

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

  ffprintf(rp.output, "\n");
  fclose(rp.output);
  delete g;
  return 0;
}
