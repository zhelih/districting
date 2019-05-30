// source file for building hess model for further usage
#include <cstdio>
#include <vector>
#include <algorithm>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"

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
  for(int i = 0; i < n; ++i)
    for(int j = 0; j < n; ++j)
      if(!F0[i][j] && !F1[i][j])
        p.h[n*i+j] = cur++; //FIXME implicit reuse of the map (i,j) -> n*i+j

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

	// add constraints (1b)
	for (int i = 0; i < n; ++i)
	{
		GRBLinExpr constr = 0;
		for (int j = 0; j < n; ++j)
			constr += X(i, j);
		model->addConstr(constr == 1);
	}

	// add constraint (1c)
	expr = 0;
	for (int j = 0; j < n; ++j)
		expr += X(j, j);
	model->addConstr(expr == k);

	// add constraint (1d)
	for (int j = 0; j < n; ++j)
	{
		GRBLinExpr constr = 0;
		for (int i = 0; i < n; ++i)
			constr += population[i] * X(i, j);
		// add for j
		model->addConstr(constr - U * X(j, j) <= 0); // U
		model->addConstr(constr - L * X(j, j) >= 0); // L
	}

	// add contraints (1e)
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (i != j)
				model->addConstr(X(i, j) <= X(j, j));

	model->update();

	//model->write("debug_hess.lp");

	return p;
}

vector<int> HessHeuristic(graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, string arg_model, double &UB, int maxIterations)
{
	vector<int> heuristicSolution(g->nr_nodes, -1);

	if (arg_model != "hess")
	{
		cerr << "Model " << arg_model << " is currently not supported by HessHeuristic.\n";
		return heuristicSolution;
	}

	vector<int> centers(k, -1);
	for (int i = 0; i < k; ++i) centers[i] = i; // dummy initialization just for build_hess_restricted. 
	vector<int> iterHeuristicSolution(g->nr_nodes, -1);

	try {

		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBVar** x = 0;
		x = build_hess_restricted(&model, g, w, population, centers, L, U, k);
		model.set(GRB_DoubleParam_TimeLimit, 60.);
		model.set(GRB_IntParam_OutputFlag, 0);

		vector<int> allNodes(g->nr_nodes, -1);
		for (int i = 0; i < g->nr_nodes; ++i) allNodes[i] = i;

		for (int iter = 0; iter < maxIterations; ++iter)
		{
			// select k centers at random
			random_shuffle(allNodes.begin(), allNodes.end());
			for (int i = 0; i < k; ++i)
				centers[i] = allNodes[i];
			
			double iterUB = INFINITY;	// the best UB found in this iteration (iter)
			double oldIterUB;			// the UB found in the previous iteration
			for (int i = 0; i < g->nr_nodes; ++i) // reset vector for this iteration's heuristic solution 
				iterHeuristicSolution[i] = -1;

			model.set(GRB_DoubleParam_MIPGap, 0.1); // Allow a loose gap in the first iteration. (The centers will be awful at first.)
			bool centersChange;

			// perform Hess descent
			do {
				oldIterUB = iterUB;
				centersChange = false;

				// Reset objective (Gurobi assumes minimization objective.)
				for (int i = 0; i < g->nr_nodes; ++i)
					for (int j = 0; j < k; ++j)
						x[i][j].set(GRB_DoubleAttr_Obj, w[i][centers[j]]);

				GRBLinExpr numCenters = 0;
				for (int i = 0; i < k; ++i)
				{
					int j = centers[i];
					numCenters += x[j][i];
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
					for (int j = 0; j < k; ++j)  // find all nodes assigned to the j-th center.
					{
						vector<int> district;
						for (int i = 0; i < g->nr_nodes; ++i)
							if (x[i][j].get(GRB_DoubleAttr_X) > 0.5)
							{
								district.push_back(i);
								iterHeuristicSolution[i] = centers[j]; // update this iteration's incumbent
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

						if (centers[j] != bestCenter) // if the centers haven't changed, there's no reason to resolve the IP
						{
							centersChange = true;
							centers[j] = bestCenter; // update the best center for this district
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
	return heuristicSolution;
}

// adds hess model constraints and the objective function to model using graph "g", distance data "dist", population data "pop"
// returns "x" variables in the Hess model
GRBVar** build_hess_restricted(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, const vector<int>&centers, int L, int U, int k)
{
	// create GUROBI Hess model
	int n = g->nr_nodes;
	int c = centers.size();
	cout << "# centers = " << c << ", while k = " << k << " and n = " << n << endl;
	if (c != k) cerr << "ERROR: improper number of centers selected for restricted problem.\n";

	// create n*(# centers) variables
	GRBVar** x = new GRBVar*[n];
	for (int i = 0; i < n; ++i)
		x[i] = model->addVars(c, GRB_BINARY);
	model->update();

	// add objective
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < c; ++j)
			x[i][j].set(GRB_DoubleAttr_Obj, w[i][centers[j]]);

	// add constraints (1b)
	for (int i = 0; i < n; ++i)
	{
		GRBLinExpr constr = 0;
		for (int j = 0; j < c; ++j)
			constr += x[i][j];
		model->addConstr(constr == 1);
	}

	//// add constraint (1c)
	cout << "Currently not adding x(V)=k constraint (fixing center vars to 1) nor coupling in build_hess_restricted." << endl;
	//GRBLinExpr expr = 0;
	//for (int j = 0; j < c; ++j)
	//	expr += x[centers[j]][j];
	//model->addConstr(expr == k);

	// add constraint (1d)
	for (int j = 0; j < c; ++j)
	{
		GRBLinExpr constr = 0;
		for (int i = 0; i < n; ++i)
			constr += population[i] * x[i][j];
		// add for j
		model->addConstr(constr - U <= 0); // U
		model->addConstr(constr - L >= 0); // L
	}

	model->update();

	return x;
}
void LocalSearch(graph* g, const vector<vector<double> >& w, const vector<int>& population,
	int L, int U, int k, vector<int>&heuristicSolution, string arg_model, double &UB) // , cvv &F0)
{
	cout << endl << "Beginning LOCAL SEARCH with UB = " << UB << "\n\n";
	if (arg_model != "hess")
	{
		cerr << "Model " << arg_model << " is currently not supported by HessHeuristic.\n";
		return;
	}
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

	//int numAvailableCenters = 0;
	//for (int i = 0; i < g->nr_nodes; ++i)
	//	if (!F0[i][i])
	//		numAvailableCenters++;

	//cout << "Number of available centers = " << numAvailableCenters << endl;
	//if (numAvailableCenters <= k) {
	//	cout << "No centers left to explore.\n";
	//	return; // no local moves possible
	//}

	try {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBVar** x = 0;
		x = build_hess_restricted(&model, g, w, population, centers, L, U, k);
		for (int i = 0; i < k; ++i)
		{
			int j = centers[i];
			x[j][i].set(GRB_DoubleAttr_LB, 1); // assign the centers to themselves.
		}
		model.set(GRB_DoubleParam_TimeLimit, 60.);
		model.set(GRB_IntParam_OutputFlag, 0);

		// create LP relaxation model. If its objective is bad, then we can terminate a local search move early. Roughly 0.2 sec (vs 2.0 sec) for MI.
		GRBModel LPmodel = GRBModel(env);
		GRBVar** LPx = 0;
		LPx = build_hess_restricted(&LPmodel, g, w, population, centers, L, U, k);
		for (int i = 0; i < g->nr_nodes; ++i)
			for (int j = 0; j < k; ++j)
			{
				LPx[i][j].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
				LPx[i][j].set(GRB_DoubleAttr_LB, 0);
				LPx[i][j].set(GRB_DoubleAttr_UB, 1);
			}
		for (int i = 0; i < k; ++i)
		{
			int j = centers[i];
			LPx[j][i].set(GRB_DoubleAttr_LB, 1); // assign the centers to themselves.
		}
		LPmodel.set(GRB_IntParam_OutputFlag, 0);

		bool improvement;
		
		do {
			improvement = false;
			for (int p = 0; p < k && !improvement; ++p)
			{
				int v = centers[p];
				cout << "  checking neighbors of node " << v << endl;

				for (int u : g->nb(v))  // swap v for u?
				{
					if (improvement) break;
					//if (F0[u][u]) continue;

					//cout << "checking node " << u << endl;

					// update cost coefficients, as if we had centers[p] = u
					for (int i = 0; i < g->nr_nodes; ++i)
						LPx[i][p].set(GRB_DoubleAttr_Obj, w[i][u]);

					LPx[v][p].set(GRB_DoubleAttr_LB, 0);
					LPx[u][p].set(GRB_DoubleAttr_LB, 1);

					LPmodel.optimize();
					double LPobj = LPmodel.get(GRB_DoubleAttr_ObjVal);

					// revert back
					for (int i = 0; i < g->nr_nodes; ++i)
						LPx[i][p].set(GRB_DoubleAttr_Obj, w[i][v]);

					LPx[v][p].set(GRB_DoubleAttr_LB, 1);
					LPx[u][p].set(GRB_DoubleAttr_LB, 0);

					//cout << "LPobj = " << LPobj << ", while UB = " << UB << endl;
					if (LPobj >= UB)
					{
						//cout << "LPobj >= UB, so terminating early.\n";
						continue; // no need to solve IP, since LP relaxation will not improve the UB.
					}

					// update cost coefficients, as if we had centers[p] = u
					for (int i = 0; i < g->nr_nodes; ++i)
						x[i][p].set(GRB_DoubleAttr_Obj, w[i][u]);

					x[v][p].set(GRB_DoubleAttr_LB, 0);
					x[u][p].set(GRB_DoubleAttr_LB, 1);

					model.optimize();

					// revert back
					for (int i = 0; i < g->nr_nodes; ++i)
						x[i][p].set(GRB_DoubleAttr_Obj, w[i][v]);

					x[v][p].set(GRB_DoubleAttr_LB, 1);
					x[u][p].set(GRB_DoubleAttr_LB, 0);

					// update incumbent (if needed)
					if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) // model was solved to optimality (subject to tolerances), so update UB.
					{
						double newUB = model.get(GRB_DoubleAttr_ObjVal);

						if (newUB < UB)
						{
							improvement = true;
							cout << "found better UB from LS restricted IP = " << newUB;
							UB = newUB;

							// update centers, costs, and var fixings
							centers[p] = u;

							for (int i = 0; i < g->nr_nodes; ++i)
								LPx[i][p].set(GRB_DoubleAttr_Obj, w[i][u]);

							LPx[v][p].set(GRB_DoubleAttr_LB, 0);
							LPx[u][p].set(GRB_DoubleAttr_LB, 1);

							for (int i = 0; i < g->nr_nodes; ++i)
								x[i][p].set(GRB_DoubleAttr_Obj, w[i][u]);

							x[v][p].set(GRB_DoubleAttr_LB, 0);
							x[u][p].set(GRB_DoubleAttr_LB, 1);

							cout << " with centers : ";
							for (int i = 0; i < k; ++i)
								cout << centers[i] << " ";
							cout << endl;

							// update heuristicSolution
							for (int i = 0; i < g->nr_nodes; ++i)
								for (int j = 0; j < k; ++j)
									if (x[i][j].get(GRB_DoubleAttr_X) > 0.5)
										heuristicSolution[i] = centers[j];
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
