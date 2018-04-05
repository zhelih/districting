#include "GRBInterface.h"
#include "KGraph.h"
#include <sstream>
#include <vector>
using namespace std;

string itos(int i) { stringstream s; s << i; return s.str(); }

string statusNumtoString(int num)
{
	if (num == 1) return "Model is loaded, but no solution information is available.";
	else if (num == 2) return "Model was solved to optimality (subject to tolerances), and an optimal solution is available.";
	else if (num == 3) return "Model was proven to be infeasible.";
	else if (num == 4) return "Model was proven to be either infeasible or unbounded. To obtain a more definitive conclusion, set the DualReductions parameter to 0 and reoptimize.";
	else if (num == 5) return "Model was proven to be unbounded. Important note: an unbounded status indicates the presence of an unbounded ray that allows the objective to improve without limit. It says nothing about whether the model has a feasible solution. If you require information on feasibility, you should set the objective to zero and reoptimize.";
	else if (num == 6) return "Optimal objective for model was proven to be worse than the value specified in the Cutoff parameter. No solution information is available.";
	else if (num == 7) return "Optimization terminated because the total number of simplex iterations performed exceeded the value specified in the IterationLimit parameter, or because the total number of barrier iterations exceeded the value specified in the BarIterLimit parameter.";
	else if (num == 8) return "Optimization terminated because the total number of branch-and-cut nodes explored exceeded the value specified in the NodeLimit parameter.";
	else if (num == 9) return "Optimization terminated because the time expended exceeded the value specified in the TimeLimit parameter.";
	else if (num == 10) return "Optimization terminated because the number of solutions found reached the value specified in the SolutionLimit parameter.";
	else if (num == 11) return "Optimization was terminated by the user.";
	else if (num == 12) return "Optimization was terminated due to unrecoverable numerical difficulties.";
	else if (num == 13) return "Unable to satisfy optimality tolerances; a sub-optimal solution is available.";
	else if (num == 14) return "An asynchronous optimization call was made, but the associated optimization run is not yet complete.";
	else if (num >= 15 || num <= 0) return "No specific error could be recognized.";
}



vector<vector<long>> solveMST(KGraph &g1, KGraph &g2)
{
	vector<vector<long>> spt;

	//bool con = g1.IsConnected();

	//cerr << con << endl;

	if (!g1.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No Spanning Tree exists! Graph is not connected." << endl;
		return spt;
	}
	try {
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		cerr << "Start Defining Y Vars" << endl;
		GRBVar **Y = new GRBVar*[g1.m];
		for (long i = 0; i < g1.m; i++)
		{
			Y[i] = model.addVars(2, GRB_BINARY);
		}
		cerr << "Start Defining Z Vars" << endl;
		GRBVar **Z = new GRBVar*[g2.m];
		for (long i = 0; i < g2.m; i++)
		{
			Z[i] = model.addVars(2, GRB_BINARY);
		}

		model.update();
		cerr << "Finish Defining Vars" << endl;

		cerr << "Adding objective function" << endl;
		for (int i = 0; i < g1.m; i++)
		{
			Y[i][0].set(GRB_DoubleAttr_Obj, 1);
			Y[i][1].set(GRB_DoubleAttr_Obj, 1);
		}

		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();

		cerr << "Adding incoming constraints in root of primal" << endl;
		long r = g1.pdadj[0][0];
		long r_d = g2.pdadj[0][0];
		cerr << "primal root: " << r << "dual root: " << r_d << endl;
		GRBLinExpr expr = 0;
		for (long j = 0; j < g1.m; j++)
		{
			if (g1.pdadj[j][0] == r && g1.pdadj[j][0] == g1.pdadj[j][1])
				expr += Y[j][0];
			else
			{
				if (g1.pdadj[j][0] == r)
					expr += Y[j][0];

				if (g1.pdadj[j][1] == r)
					expr += Y[j][1];
			}
		}
		model.addConstr(expr == 0);
		model.update();
		cerr << "Adding incoming constraints in non-root of primal" << endl;

		cerr << "nodes of primal: " << g1.n << endl;

		for (long i = 0; i < g1.n; i++)
		{
			if (i == r)
				continue;
			GRBLinExpr expr1 = 0;
			for (long j = 0; j < g2.m; j++)
			{
				if (g1.pdadj[j][0] == i && g1.pdadj[j][0] == g1.pdadj[j][1])
					expr1 += Y[j][0];
				else
				{
					if (g1.pdadj[j][0] == i)
						expr1 += Y[j][0];

					else if (g1.pdadj[j][1] == i)
						expr1 += Y[j][1];
				}

			}
			model.addConstr(expr1 == 1);
		}
		model.update();

		cerr << "Adding incoming constraints in root of dual" << endl;
		GRBLinExpr expr2 = 0;
		for (long j = 0; j < g2.m; j++)
		{
			if (g2.pdadj[j][0] == r_d && g2.pdadj[j][0] == g2.pdadj[j][1])
				expr2 += Z[j][0];

			else
			{
				if (g2.pdadj[j][0] == r_d)
					expr2 += Z[j][0];
				if (g2.pdadj[j][1] == r_d)
					expr2 += Z[j][1];
			}
		}
		model.addConstr(expr2 == 0);
		model.update();
		cerr << "Adding incoming constraints in non-root of dual" << endl;

		for (long i = 0; i < g2.n; i++)
		{
			if (i == r_d)
				continue;
			GRBLinExpr expr3 = 0;
			for (long j = 0; j < g2.m; j++)
			{
				if (g2.pdadj[j][0] == i && g2.pdadj[j][0] == g2.pdadj[j][1])
					expr3 += Z[j][0];
				else
				{
					if (g2.pdadj[j][0] == i)
						expr3 += Z[j][0];
					if (g2.pdadj[j][1] == i)
						expr3 += Z[j][1];
				}

			}
			model.addConstr(expr3 == 1);
		}
		model.update();

		cerr << "Adding intersection constraints" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			GRBLinExpr expr4 = 0;
			if (g1.pdadj[i][0] == g1.pdadj[i][1])
			{
				expr4 = Y[i][0] + Z[i][0] + Z[i][1];
			}
			else if (g2.pdadj[i][0] == g2.pdadj[i][1])
			{
				expr4 = Y[i][0] + Y[i][1] + Z[i][0];
			}
			else
			{
				expr4 = Y[i][0] + Y[i][1] + Z[i][0] + Z[i][1];
			}
			model.addConstr(expr4 == 1);
		}
		model.update();
		model.write("debug.lp");
		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;


		for (int i = 0; i < g1.m; i++)
		{
			vector<long> art;
			if (Y[i][0].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				art.push_back(g1.pdadj[i][0]);
			}
			if (Y[i][1].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				art.push_back(g1.pdadj[i][1]);
			}
			spt.push_back(art);
		}

		delete[] Y;
		delete[] Z;


	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return spt;
}
