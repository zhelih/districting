// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include "LazyConstraints.h"
#include <chrono>

void build_cut1(GRBModel* model, GRBVar** x, graph* g)
{
	// create n^2 variables, and set UB=0
	int n = g->nr_nodes;
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
	GRBVar** y = new GRBVar*[n]; //FIXME ever deleted?
	for (int i = 0; i < n; ++i)
	{
		GRBVar *y_temp = new GRBVar[n];
		for (int j = 0; j < n; ++j)
		{
			// set UBs to zero for all y vars
			y_temp[j] = model->addVar(0.0, 0.0, 0.0, GRB_BINARY);
		}
		y[i] = y_temp;
	}

	// add constraint (23b)
	for (int i = 0; i < n; ++i)
	{
		GRBLinExpr expr = 0;
		for (int j : g->nb(i))
		{
			// set UBs to one for all arcs
			y[j][i].set(GRB_DoubleAttr_UB, 1.0);
			expr += y[j][i];
		}
		model->addConstr(expr + x[i][i] == 1);
	}
  //FIXME delete
	LazyConstraints* cb = new LazyConstraints(x, y, g);	// tell Gurobi which function generates the lazy cuts.
	model->setCallback(cb);
	model->update();
	//optimize the model
	//model->optimize();

	cerr << endl;
	cerr << "Number of callbacks : " << LazyConstraints::numCallbacks << endl;
	cerr << "Time in callbacks : " << LazyConstraints::TotalCallbackTime << endl;
	cerr << "Number of lazy cuts : " << LazyConstraints::numLazyCuts << endl;
}

void Cut2Callback::callback()
{
	using namespace std;
	try
	{
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;
			auto start = chrono::steady_clock::now();

			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					x[i][j] = getSolution(grb_x[i][j]);

			// try clusterheads
			bool done = false;
			for (int i = 0; i < n && !done; ++i)
			{
				if (x[i][i] > 0.5) // i is a clusterhead
				{
					// run DFS from i on C_i, compute A(C_i) to save time later

					// clean A(C_i) and visited
					for (int j = 0; j < n; ++j)
					{
						aci[j] = 0;
						visited[j] = 0;
					}

					// run DFS on C_i
					s.clear(); s.push_back(i);
					while (!s.empty())
					{
						int cur = s.back(); s.pop_back();
						visited[cur] = 1;
						for (int nb_cur : g->nb(cur))
							if (x[nb_cur][i] > 0.5) // if nb_cur is in C_i
							{
								if (!visited[nb_cur])
									s.push_back(nb_cur);
							}
							else aci[nb_cur] = 1; // nb_cur is a neighbor of a vertex in C_i, thus in A(C_i)a
					}

					// here if C_i is connected, all vertices in C_i must be visited
					for (int j = 0; j < n; ++j)
						if (x[j][i] > 0.5 && !visited[j])
						{
							// compute i-j separator, A(C_i) is already computed)
							GRBLinExpr expr = 0;
							// start DFS from j to find R_i
							for (int k = 0; k < n; ++k)
								visited[k] = 0;
							s.clear(); s.push_back(j); visited[j] = 1;
							while (!s.empty())
							{
								int cur = s.back(); s.pop_back();
								for (int nb_cur : g->nb(cur))
								{
									if (!visited[nb_cur])
									{
										visited[nb_cur] = 1;
										if (aci[nb_cur])
											expr += grb_x[nb_cur][i];
										else s.push_back(nb_cur);
									}
								}
							}
							expr -= grb_x[j][i]; // RHS
							addLazy(expr >= 0);
							done = true;
							break;
						}
				}
			}
			chrono::duration<double> d = chrono::steady_clock::now() - start;
			callbackTime += d.count();
		}
	}
	catch (GRBException e)
	{
		fprintf(stderr, "Error number: %d\n", e.getErrorCode());
		fprintf(stderr, "%s\n", e.getMessage().c_str());
	}
	catch (...)
	{
		fprintf(stderr, "Error during callback\n");
	}
}

Cut2Callback* build_cut2(GRBModel* model, GRBVar** x, graph* g)
{
	model->getEnv().set(GRB_IntParam_LazyConstraints, 1); // turns off presolve!!!
	Cut2Callback* cb = new Cut2Callback(x, g);
	model->setCallback(cb);
	return cb;
}
