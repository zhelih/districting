// source file for building hess model for further usage
#include <cstdio>
#include <vector>
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
GRBVar** build_hess(GRBModel* model, graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, cvv& F0, cvv& F1)
{
	// create GUROBI Hess model
	int n = g->nr_nodes;

	// create n^2 variables
	GRBVar** x = new GRBVar*[n];
	for (int i = 0; i < n; ++i)
	{
		x[i] = new GRBVar[n];
		for (int j = 0; j < n; ++j)
			if (IS_X(i, j))
				x[i][j] = model->addVar(0., 1., 0., GRB_BINARY);
	}
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

	return x;
}

vector<int> HessHeuristic(graph* g, const vector<vector<double> >& w, const vector<int>& population, int L, int U, int k, vector<int>&centers, string arg_model, double &UB)
{
	if (centers.empty())
	{
		cout << "Set of centers is currently empty. Using the set { 0, 1, ..., k-1 }.\n";
		centers.resize(k, 0);
		for (int i = 0; i < k; ++i)
			centers[i] = i;
	}

	vector<int> heuristicSolution(g->nr_nodes, 0);

	if (arg_model != "hess")
	{
		cerr << "Model " << arg_model << " is currently not supported by HessHeuristic.\n";
		return heuristicSolution;
	}

	double newUB = INFINITY;
	bool CentersChanged;

	try {

		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		GRBVar** x = 0;
		x = build_hess_restricted(&model, g, w, population, centers, L, U, k);
		model.set(GRB_DoubleParam_TimeLimit, 600.); 

		do {
			UB = newUB;  // the UBs are monotone, so this is okay to do.
			CentersChanged = false;

			// Set objective (Gurobi assumes minimization objective.)
			for (int i = 0; i < g->nr_nodes; ++i)
				for (int j = 0; j < k; ++j)
					x[i][j].set(GRB_DoubleAttr_Obj, w[i][centers[j]]);

			model.optimize();

			if (model.get(GRB_IntAttr_Status) == 2 || model.get(GRB_IntAttr_Status) == 9) // model was solved to optimality (subject to tolerances), so update UB.
			{
				newUB = model.get(GRB_DoubleAttr_ObjVal);
				cout << "UB from restricted IP = " << newUB << endl;
			}
			if (newUB < UB) // found a better solution; update the centers for the next iteration.
			{
				for (int j = 0; j < k; ++j)
				{
					vector<int> district;
					for (int i = 0; i < g->nr_nodes; ++i)
						if (x[i][j].get(GRB_DoubleAttr_X) > 0.5)
						{
							district.push_back(i);
							heuristicSolution[i] = centers[j];
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
							//Cost += ((double)dist[v][c] / 1000.) * ((double)dist[v][c] / 1000.) * population[v];
						}

						if (Cost < bestCost)
						{
							bestCenter = c;
							bestCost = Cost;
						}
					}

					if (centers[j] != bestCenter)
					{
						CentersChanged = true;
						centers[j] = bestCenter; // update the best center for this district
					}
				}
			}
		} while (newUB < UB && CentersChanged); // if the centers did not change, then the next iteration will not improve the UB, so terminate early.

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
	
	UB = newUB;
	cerr << "UB at end of HessHeuristic = " << UB << endl;

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

	// add constraints (1b)
	for (int i = 0; i < n; ++i)
	{
		GRBLinExpr constr = 0;
		for (int j = 0; j < c; ++j)
			constr += x[i][j];
		model->addConstr(constr == 1);
	}

	// add constraint (1c)
	GRBLinExpr expr = 0;
	for (int j = 0; j < c; ++j)
		expr += x[centers[j]][j];
	model->addConstr(expr == k);

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