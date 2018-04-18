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
		long r = g1.edge[0][0];
		long r_d = g2.edge[0][0];
		long u, v;
		cerr << "primal root: " << r << "dual root: " << r_d << endl;
		vector <GRBLinExpr> exprPrim(g1.n,0);
		vector <GRBLinExpr> exprDual(g2.n, 0);

		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			exprPrim[u] += Y[i][0];
			exprPrim[v] += Y[i][1];
			u = g2.edge[i][0];
			v = g2.edge[i][1];
			exprDual[u] += Z[i][0];
			exprDual[v] += Z[i][1];
		}

		cerr << "Adding incoming constraints of primal" << endl;

		model.addConstr(exprPrim[r] == 0);
		
		for (long i = 0; i < g1.n; i++)
		{
			if (i == r)
				continue;
			model.addConstr(exprPrim[i] == 1);
		}
		model.update();

		cerr << "Adding incoming constraints of dual" << endl;
		model.addConstr(exprDual[r_d] == 0);

		for (long i = 0; i < g2.n; i++)
		{
			if (i == r_d)
				continue;
			model.addConstr(exprDual[i] == 1);
		}
		model.update();

		cerr << "Adding intersection constraints" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			model.addConstr(Y[i][0] + Y[i][1] + Z[i][0] + Z[i][1] == 1);
		}
		model.update();

		cerr << "Trimming loops" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			if (u == v)
				Y[i][1].set(GRB_DoubleAttr_UB, 0.0);
			u = g2.edge[i][0];
			v = g2.edge[i][1];
			if (u == v)
				Z[i][1].set(GRB_DoubleAttr_UB, 0.0);
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
		vector<long> arc(2,0);
		for (int i = 0; i < g1.m; i++)
		{
			if (Y[i][0].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				arc[0]=i;
				arc[1] = 0;
				spt.push_back(arc);
			}
			if (Y[i][1].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				arc[0] = i;
				arc[1] = 1;
				spt.push_back(arc);
			}		
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

vector<vector<long>> solveFM2(KGraph &g1, KGraph &g2, KGraph &d, long &K, long &L, long &U)
{
	vector<vector<long>> form2;
	if (!g1.IsConnected()) // check if problem is feasible before sending to Gurobi
	{
		cerr << "No Spanning Tree exists! Graph is not connected." << endl;
		return form2;
	}
	try {
		GRBEnv env = GRBEnv();
		//env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		GRBModel model = GRBModel(env);
		env.set(GRB_IntParam_Method, 3);
		//model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		cerr << "Start Defining Y, Z, X, and F Vars" << endl; //For Primal Spanning Tree
		GRBVar **Y = new GRBVar*[g1.m];
		GRBVar **Z = new GRBVar*[g2.m];
		GRBVar **X = new GRBVar*[g1.m];
		GRBVar **F = new GRBVar*[g1.m];
		for (long i = 0; i < g1.m; i++)
		{
			Y[i] = model.addVars(2, GRB_CONTINUOUS);
			Z[i] = model.addVars(2, GRB_CONTINUOUS);
			X[i] = model.addVars(2, GRB_BINARY);
			F[i] = model.addVars(2, GRB_CONTINUOUS);
		}

		cerr << "Start Defining R Vars" << endl; //For Root Variables
		GRBVar *R = model.addVars(g1.n, GRB_BINARY);

		cerr << "Start Defining G Vars" << endl; //For Generation Variables
		GRBVar *G = model.addVars(g1.n, GRB_CONTINUOUS);

		cerr << "Increase Branch Priority of R Vars" << endl; //For Root Var Priority

		for (long i = 0; i < g1.n; i++)
		{
			R[i].set(GRB_IntAttr_BranchPriority, 9.0);
		}

		model.update();
		cerr << "Finish Defining Vars" << endl;

		cerr << "Currently use hop-based distances" << endl;
		vector <vector<long>> dis (g1.m, vector<long>(2));
		for (long i = 0; i < g1.m; i++)
		{
			dis[i][0] = 1;
			dis[i][1] = 1;		
		}

		cerr << "Adding objective function" << endl;
		for (int i = 0; i < g1.m; i++)
		{
			F[i][0].set(GRB_DoubleAttr_Obj, dis[i][0]);
			F[i][1].set(GRB_DoubleAttr_Obj, dis[i][1]);
		}

		model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
		model.update();
		long r = g1.edge[0][0];
		long r_d = g2.edge[0][0];
		long u, v;
		cerr << "primal root: " << r << "dual root: " << r_d << endl;
		vector <GRBLinExpr> exprPrim(g1.n, 0);
		vector <GRBLinExpr> exprDual(g2.n, 0);

		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			exprPrim[u] += Y[i][0];
			exprPrim[v] += Y[i][1];
			u = g2.edge[i][0];
			v = g2.edge[i][1];
			exprDual[u] += Z[i][0];
			exprDual[v] += Z[i][1];
		}

		cerr << "Adding incoming constraints of primal" << endl;

		model.addConstr(exprPrim[r] == 0);

		for (long i = 0; i < g1.n; i++)
		{
			if (i == r)
				continue;
			model.addConstr(exprPrim[i] == 1);
		}
		model.update();

		cerr << "Adding incoming constraints of dual" << endl;
		model.addConstr(exprDual[r_d] == 0);

		for (long i = 0; i < g2.n; i++)
		{
			if (i == r_d)
				continue;
			model.addConstr(exprDual[i] == 1);
		}
		model.update();

		cerr << "Adding intersection constraints" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			model.addConstr(Y[i][0] + Y[i][1] + Z[i][0] + Z[i][1] == 1);
		}
		model.update();

		cerr << "Adding X <= W constraints" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			model.addConstr(Y[i][0] + Y[i][1] >= X[i][0] + X[i][1]);
		}
		model.update();


		cerr << "Adding Incoming + Root constraints" << endl; 
		vector <GRBLinExpr> exprPrimRoot(g1.n,0);

		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			exprPrimRoot[u] += X[i][0];
			exprPrimRoot[v] += X[i][1];
		}
		for (long i = 0; i < g1.n; i++)
		{
			model.addConstr(exprPrimRoot[i] + R[i] == 1);
		}
		model.update();

		

		cerr << "Adding r(V) = k constraints" << endl;
		GRBLinExpr expr6 = 0;
		for (long i = 0; i < g1.n; i++)
		{
			expr6 += R[i];
		}
		model.addConstr(expr6 == K);
		model.update();

		cerr << "Adding generation bounds' constraints" << endl;
		for (long i = 0; i < g1.n; i++)
		{
			model.addConstr(G[i] <= U*R[i]);
			model.addConstr(G[i] >= L*R[i]);
		}
		model.update();

		cerr << "Adding conservation's constraints" << endl;
		vector <GRBLinExpr> exprIncom(g1.n, 0);
		vector <GRBLinExpr> exprOutgo(g1.n, 0);

		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			exprIncom[u] += F[i][0];
			exprIncom[v] += F[i][1];
			exprOutgo[u] += F[i][1];
			exprOutgo[v] += F[i][0];
		}

		for (long i = 0; i < g1.n; i++)
		{
			model.addConstr(exprIncom[i] + G[i] == exprOutgo[i] + d.pop[i]);
		}
		model.update();

		cerr << "Adding flow bounds' constraints" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			model.addConstr(F[i][0] <= (U - d.pop[v])*X[i][0]);
			model.addConstr(F[i][0] >= d.pop[u] * X[i][0]);
			model.addConstr(F[i][1] <= (U - d.pop[u])*X[i][1]);
			model.addConstr(F[i][1] >= d.pop[v] * X[i][1]);
		}
		model.update();

		cerr << "Trimming loops" << endl;
		for (long i = 0; i < g1.m; i++)
		{
			u = g1.edge[i][0];
			v = g1.edge[i][1];
			if (u == v)
			{
				Y[i][1].set(GRB_DoubleAttr_UB, 0.0);
				X[i][1].set(GRB_DoubleAttr_UB, 0.0);
			}				
			u = g2.edge[i][0];
			v = g2.edge[i][1];
			if (u == v)
				Z[i][1].set(GRB_DoubleAttr_UB, 0.0);
		}
		model.update();
		model.write("debug.lp");

		//cerr << "Obtaining a feasible upper bound" << endl;

		//vector<long> sortedPop(g1.n,0);

		//for (long i = 0; i < g1.n; i++)
		//{
		//	sortedPop[i] = d.pop[i];
		//}
		//sort(sortedPop.begin(), sortedPop.end());
		////long u, v;
		//for (long i = g1.n - 1; i > g1.n - K -1; i--)
		//{
		//	u = sortedPop[i];
		//	for (long j = 0; j < g1.n; j++)
		//	{
		//		if (u == d.pop[j])
		//			R[j].set(GRB_DoubleAttr_Start, 1.0);
		//	}
		//}
		//model.update();

		cerr << "Optimizing" << endl;
		model.optimize();

		long NumBBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		double lb = model.get(GRB_DoubleAttr_ObjBound);
		cerr << "MIP LB = " << lb << endl;
		cerr << "Number of branch-and-bound nodes : " << NumBBNodes << endl;

		cerr << "Extracting solution" << endl;
		int status = model.get(GRB_IntAttr_Status);
		cerr << statusNumtoString(status) << endl;
		vector<long> arc (2);

		for (int i = 0; i < g1.m; i++)
		{
			
			if (X[i][0].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				arc[0] = i;
				arc[1] = 0;
				form2.push_back(arc);
			}
			if (X[i][1].get(GRB_DoubleAttr_X) > 0.5)	// since Gurobi uses floating point numbers, a binary variable may take value 0.9999. So just check if it's above 0.5
			{
				arc[0] = i;
				arc[1] = 1;
				form2.push_back(arc);
			}
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
	return form2;
}
