// source file for building hess model for further usage
#include <cstdio>
#include <vector>
#include <algorithm>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include "io.h"

using namespace std;

double get_objective_coefficient(const vector<vector<int>>& dist, const vector<int>& population, int i, int j)
{
  return (static_cast<double>(dist[i][j]) / 1000.) * (static_cast<double>(dist[i][j]) / 1000.) * static_cast<double>(population[i]);
}


// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
hess_params build_hess(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1)
{
  // create GUROBI Hess model
  int n = g->nr_nodes;
  hess_params p;
  p.n = n;
  p.F0 = F0; p.F1 = F1; // copy

  // hash variables
  int cur = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (!F0[i][j] && !F1[i][j])
        p.h[n*i + j] = cur++; //FIXME implicit reuse of the map (i,j) -> n*i+j

  printf("Build hess : created %lu variables\n", p.h.size());
  int nr_var = static_cast<int>(p.h.size());

  // create variables
  p.x = model->addVars(nr_var, GRB_BINARY);
  model->update();

  // Set objective: minimize sum d^2_ij*x_ij
  GRBLinExpr expr = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      expr += w[i][j] * X(i, j);

  model->setObjective(expr, GRB_MINIMIZE);

  // add constraints (b)
  for (int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for (int j = 0; j < n; ++j)
      constr += X(i, j);
    model->addConstr(constr == 1);
  }

  // add constraint (c)
  expr = 0;
  for (int j = 0; j < n; ++j)
    expr += X(j, j);
  model->addConstr(expr == k);

  // add constraint (d)
  for (int j = 0; j < n; ++j)
  {
    GRBLinExpr constr = 0;
    for (int i = 0; i < n; ++i)
      constr += population[i] * X(i, j);
    // add for j
    model->addConstr(constr - U * X(j, j) <= 0); // U
    model->addConstr(constr - L * X(j, j) >= 0); // L
  }

  // add contraints (e)
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      if (i != j && !F0[i][j])
        model->addConstr(X(i, j) <= X(j, j));

  model->update();

  //model->write("debug_hess.lp");

  return p;
}

// populate F0 and F1 depending on current centers
void populate_hess_params(hess_params& p, graph* g, const vector<int>& centers)
{
  int n = g->nr_nodes; p.n = n;
  int cur = 0;
  p.h.clear();

  //define x_ij for every for every j \in centers
  for (int j : centers)
    for (int i = 0; i < n; ++i)
      p.h[n*i + j] = cur++;

  // clear F0 and F1
  p.F0.resize(n); p.F1.resize(n);
  for (int i = 0; i < n; ++i)
  {
    p.F0[i].resize(n); p.F1[i].resize(n);
    fill(p.F0[i].begin(), p.F0[i].end(), false);
    fill(p.F1[i].begin(), p.F1[i].end(), false);
  }

  // fill F0 and F1
  //for(int j : centers)
  //  p.F1[j][j] = true; // the centers fixed to 1
  vector<bool> aux_F1(n, false);
  for (int j : centers)
    aux_F1[j] = true;

  for (int j = 0; j < n; ++j)
    if (!aux_F1[j])
      for (int i = 0; i < n; ++i)
        p.F0[i][j] = true; // other centers as well as corresponding i fixed to 0
}

#define ENSURE(i,j) {if(p.h.count(p.n*i+j)==0){fprintf(stderr,"ensure failed at line %d for i = %d, j = %d\n", __LINE__, i, j);exit(1);}}

// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
hess_params build_hess_restricted(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, const vector<int>&centers, int L, int U, int k)
{
  // create GUROBI Hess model
  int n = g->nr_nodes;
  int c = centers.size();
  cout << "# centers = " << c << ", while k = " << k << " and n = " << n << endl;
  if (c != k) cerr << "ERROR: improper number of centers selected for restricted problem." << endl;

  // hash (i,j) for j in centers
  hess_params p;
  populate_hess_params(p, g, centers);
  // create (n-1)*(# centers) variables
  p.x = model->addVars(p.h.size(), GRB_BINARY); // h.size() = n * centers.size()
  model->update();

  // recompute objective
  for (int i = 0; i < n; ++i)
    for (int j : centers)
    {
      ENSURE(i, j);
      X_V(i, j).set(GRB_DoubleAttr_Obj, w[i][j]);
    }

  // add constraints (1b)
  for (int i = 0; i < n; ++i)
  {
    GRBLinExpr constr = 0;
    for (int j : centers)
    {
      ENSURE(i, j);
      constr += X(i, j);
    }
    model->addConstr(constr == 1); // each i must be assigned to a center
  }

  // add constraint (1d)
  for (int j : centers)
  {
    GRBLinExpr constr = 0;
    for (int i = 0; i < n; ++i)
    {
      ENSURE(i, j);
      constr += population[i] * X(i, j);
    }
    // add for j
    model->addConstr(constr - U <= 0); // U
    model->addConstr(constr - L >= 0); // L
  }

  model->update();

  return p;
}

void ContiguityHeuristic(vector<int> &heuristicSolution, graph* g, const vector<vector<double> > &w,
  const vector<int> &population, int L, int U, int k, double &UB, string arg_model)
{
  vector<int> centers(k); 
  int pos = 0;
  for (int i = 0; i < g->nr_nodes; ++i)
  {
    if (heuristicSolution[i] == i)
    {
      centers[pos] = i;
      pos++;
      cout << "i=" << i << endl;
    }
  }

  HessCallback* cb = nullptr;

  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 600.);
    //model.set(GRB_IntParam_OutputFlag, 0);

    hess_params p = build_hess_restricted(&model, g, w, population, centers, L, U, k);
    populate_hess_params(p, g, centers); // now just use X(i,j), F0, F1 and more importantly hashtable are properly set up


    // begin zero-pop cuts
    //if (arg_model == "cut" || arg_model == "lcut")
    //{
      /*for (int i = 0; i < g->nr_nodes; ++i)
      {
        if (population[i] > 0) continue;
        cout << "population of node " << i << " is zero " << endl;
        for (int j : centers)
        {
          GRBLinExpr expr = 0;
          for (int nb : g->nb(i))
            expr += X_V(nb, j);
          model.addConstr(X_V(i, j) <= expr);
        }
      }*/
      //model.set(GRB_IntParam_Cuts, 3); // set aggressive cuts
    //}
    // end zero-pop cuts

    if (arg_model == "shir")
      build_shir(&model, p, g);
    else if (arg_model == "mcf")
      build_mcf(&model, p, g);
    else if (arg_model == "cut")
      cb = build_cut(&model, p, g, population);
    else if (arg_model == "lcut")
      cb = build_lcut(&model, p, g, population, U);
    else {
      fprintf(stderr, "ERROR: Unknown contiguity model : %s\n", arg_model.c_str());
      exit(1);
    }
    for (int i = 0; i < k; ++i)
    {
      int v = centers[i];
      X_V(v,v).set(GRB_DoubleAttr_LB, 1);
    }
    model.optimize();


    if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) // model was solved to optimality (subject to tolerances), so update UB.
    {
      UB = model.get(GRB_DoubleAttr_ObjVal); // we get a warm start from previous round, so iterUB will only get better
      cout << "  UB from ContiguityHeuristic restricted IP = " << UB << " using centers : ";
      for (int i = 0; i < k; ++i)
        cout << centers[i] << " ";
      cout << endl;

      for (int i = 0; i < g->nr_nodes; ++i)
        for (int j = 0; j < g->nr_nodes; ++j)
          if (p.F0[i][j]) continue; else if (p.F1[i][j] || X_V(i, j).get(GRB_DoubleAttr_X) > 0.5)
            heuristicSolution[i] = j;
    }
  }
  catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  }
  catch (const char* msg) {
    cout << "Exception with message : " << msg << endl;
  }
  catch (...) {
    cout << "Exception during optimization" << endl;
  }

  cout << "UB at end of ContiguityHeuristic = " << UB << endl;
  double obj = 0;
  for (int i = 0; i < g->nr_nodes; ++i)
    obj += w[i][heuristicSolution[i]];
  cout << "UB of (contiguous) heuristicSolution = " << obj << endl;
  if (cb)
    delete cb;
  return;
}

vector<int> HessHeuristic(graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, double &UB, int maxIterations, bool do_cuts)
{
  vector<int> heuristicSolution(g->nr_nodes, -1);

  vector<int> centers(k, -1);
  for (int i = 0; i < k; ++i) centers[i] = i; // dummy initialization just for build hess restricted. 
  vector<int> iterHeuristicSolution(g->nr_nodes, -1);

  HessCallback* cb = nullptr;

  try {

    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 60.);
    model.set(GRB_IntParam_OutputFlag, 0);

    hess_params p = build_hess_restricted(&model, g, w, population, centers, L, U, k);

    vector<int> allNodes(g->nr_nodes, -1);
    for (int i = 0; i < g->nr_nodes; ++i) allNodes[i] = i;

    for (int iter = 0; iter < maxIterations; ++iter)
    {
      // select k centers at random
      random_shuffle(allNodes.begin(), allNodes.end()); //FIXME O(n) instead of O(k)
      for (int i = 0; i < k; ++i)
        centers[i] = allNodes[i];

      double iterUB = INFINITY; // the best UB found in this iteration (iter)
      double oldIterUB;     // the UB found in the previous iteration
      for (int i = 0; i < g->nr_nodes; ++i) // reset vector for this iteration's heuristic solution 
        iterHeuristicSolution[i] = -1;

      model.set(GRB_DoubleParam_MIPGap, 0.1); // Allow a loose gap in the first iteration. (The centers will be awful at first.)
      bool centersChange;

      // perform Hess descent
      do {
        oldIterUB = iterUB;
        centersChange = false;

        // Reset objective (Gurobi assumes minimization objective.)
        populate_hess_params(p, g, centers); // now just use X(i,j), F0, F1 and more importantly hashtable are properly set up
        if (do_cuts)
        {
          if (cb)
            delete cb;
          build_cut(&model, p, g, population); //FIXME do pointers instead? worth it? prob no
        }
        model.reset(); // should be done in any case for predicted behavior
        for (int i = 0; i < g->nr_nodes; ++i)
          for (int j : centers)
          {
            ENSURE(i, j);
            X_V(i, j).set(GRB_DoubleAttr_Obj, w[i][j]);
          }

        GRBLinExpr numCenters = 0;
        for (int j : centers)
        {
          ENSURE(j, j);
          numCenters += X(j, j); // different map for different centers!
        }

        model.addConstr(numCenters == k, "fixCenters");  // in essence, fix the current centers to 1
        model.optimize();

        if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) // model was solved to optimality (subject to tolerances), so update UB.
        {
          iterUB = model.get(GRB_DoubleAttr_ObjVal); // we get a warm start from previous round, so iterUB will only get better
          cout << "  UB from restricted IP = " << iterUB << " using centers : ";
          for (int i = 0; i < k; ++i)
            cout << centers[i] << " ";
          cout << endl;
        }

        model.remove(model.getConstrByName("fixCenters")); // unfix the current centers

        if (iterUB < oldIterUB) // objective value strictly improved, so we need to update incumbent for this iteration
        {
          for (int j_i = 0; j_i < k; ++j_i)  // find all nodes assigned to the j-th center.
          {
            vector<int> district;
            for (int i = 0; i < g->nr_nodes; ++i)
            {
              ENSURE(i, centers[j_i]);
              if (X_V(i, centers[j_i]).get(GRB_DoubleAttr_X) > 0.5) //FIXME cerr if not IS_X
              {
                district.push_back(i);
                iterHeuristicSolution[i] = centers[j_i]; // update this iteration's incumbent
              }
            }

            // find best center of this district
            int bestCenter = -1;
            double bestCost = INFINITY;

            for (int i = 0; i < district.size(); ++i)
            {
              // compute the Cost of this district when centered at c
              int c = district[i];
              double Cost = 0;

              for (int p = 0; p < district.size(); ++p)
              {
                int v = district[p];
                Cost += w[v][c];
              }

              if (Cost < bestCost)
              {
                bestCenter = c;
                bestCost = Cost;
              }
            }

            if (centers[j_i] != bestCenter) // if the centers haven't changed, there's no reason to resolve the IP
            {
              centersChange = true;
              centers[j_i] = bestCenter; // update the best center for this district
            }
          }
        }
        // Now that the centers should be reasonable, require a tighter tolerance
        model.set(GRB_DoubleParam_MIPGap, 0.0005);

      } while (iterUB < oldIterUB && centersChange);

      // update incumbents (if needed)
      if (iterUB < UB)
      {
        UB = iterUB;
        heuristicSolution = iterHeuristicSolution;
      }
      cout << "In iteration " << iter << " of HessHeuristic, objective value of incumbent is = " << UB << endl;
    }
  }
  catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  }
  catch (const char* msg) {
    cout << "Exception with message : " << msg << endl;
  }
  catch (...) {
    cout << "Exception during optimization" << endl;
  }

  cout << "UB at end of HessHeuristic = " << UB << endl;
  double obj = 0;
  for (int i = 0; i < g->nr_nodes; ++i)
    obj += w[i][heuristicSolution[i]];
  cout << "UB of heuristicSolution = " << obj << endl;
  if (cb)
    delete cb;
  return heuristicSolution;
}

void LocalSearch(graph* g, const vector<vector<double> >& w, const vector<int>& population,
  int L, int U, int k, vector<int>&heuristicSolution, string arg_model, double &UB) // , cvv &F0)
{
    cout << endl << "Beginning LOCAL SEARCH with UB = " << UB << "\n\n";
    // initialize the centers from heuristicSolution
    vector<int> centers(k, -1);
    int pos = 0;
    for (int i = 0; i < g->nr_nodes; ++i)
    {
      if (heuristicSolution[i] == i)
      {
        centers[pos] = i;
        pos++;
      }
    }
    try {
      GRBEnv env = GRBEnv();
      GRBModel model = GRBModel(env);
      hess_params p = build_hess_restricted(&model, g, w, population, centers, L, U, k);
      for (int j : centers)
      {
        ENSURE(j, j);
        X_V(j, j).set(GRB_DoubleAttr_LB, 1); // assign the centers to themselves.
      }
      model.set(GRB_DoubleParam_TimeLimit, 60.);
      model.set(GRB_IntParam_OutputFlag, 0);
      bool improvement;
      do {
        improvement = false;
        for (int c_i = 0; c_i < k && !improvement; ++c_i)
        {
          int v = centers[c_i];
          cout << "  checking neighbors of node " << v << endl;
          for (int u : g->nb(v))  // swap v for u?
          {
            if (improvement) break;
            model.reset();
            model.set(GRB_DoubleParam_Cutoff, UB);
            // update cost coefficients, as if we had centers[p] = u
            for (int i = 0; i < g->nr_nodes; ++i)
            {
              ENSURE(i,v);
              X_V(i,v).set(GRB_DoubleAttr_Obj, w[i][u]);
            }
            X_V(v,v).set(GRB_DoubleAttr_LB, 0);
            X_V(u,v).set(GRB_DoubleAttr_LB, 1);
            model.optimize();
            // revert back
            for (int i = 0; i < g->nr_nodes; ++i)
              X_V(i,v).set(GRB_DoubleAttr_Obj, w[i][v]);
            X_V(v,v).set(GRB_DoubleAttr_LB, 1);
            X_V(u,v).set(GRB_DoubleAttr_LB, 0);
            // update incumbent (if needed) if solved or timed out
            if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) // model was solved to optimality (subject to tolerances), so update UB.
            {
              double newUB = model.get(GRB_DoubleAttr_ObjVal);
              if (newUB < UB)
              {
                improvement = true;
                cout << "found better UB from LS restricted IP = " << newUB;
                UB = newUB;
                // update centers, costs, and var fixings
                centers[c_i] = u;
                for (int i = 0; i < g->nr_nodes; ++i)
                  X_V(i,v).set(GRB_DoubleAttr_Obj, w[i][u]);
                X_V(v,v).set(GRB_DoubleAttr_LB, 0);
                X_V(u,v).set(GRB_DoubleAttr_LB, 1);
                cout << " with centers : ";
                for (int i = 0; i < k; ++i)
                  cout << centers[i] << " ";
                cout << endl;
                // update hess params
                populate_hess_params(p, g, centers);
                // update heuristicSolution
                for (int i = 0; i < g->nr_nodes; ++i)
                  for (int j = 0; j < k; ++j)
                  {
                    ENSURE(i,centers[j]);
                    if (X_V(i,centers[j]).get(GRB_DoubleAttr_X) > 0.5)
                      heuristicSolution[i] = centers[j];
                  }
              }
            }
          }
        }
      } while (improvement);
    }
    catch (GRBException e) {
      cout << "Error code = " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    }
    catch (const char* msg) {
      cout << "Exception with message : " << msg << endl;
    }
    catch (...) {
      cout << "Exception during optimization" << endl;
    }
    cout << "UB at end of local search heuristic = " << UB << endl;
    double obj = 0;
    for (int i = 0; i < g->nr_nodes; ++i)
      obj += w[i][heuristicSolution[i]];
    cout << "UB of heuristicSolution = " << obj << endl;
}
