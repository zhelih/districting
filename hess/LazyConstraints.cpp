#include "LazyConstraints.h"
#include "gurobi_c++.h"
#include "graph.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <queue>
#include <set>
#include <map>
#include <omp.h>
#include <memory>
#include <sstream>
#include <vector>
using namespace std;


LazyConstraints::LazyConstraints(GRBVar **xvars, GRBVar **yvars, graph *g): g1(g)
{
	varsx = xvars;
	varsy = yvars;
  numCallbacks = 0L;
  TotalCallbackTime = 0.;
  numLazyCuts = 0L;
};

void LazyConstraints::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			numCallbacks++;
			time_t start = clock();

			double **X = new double*[g1->nr_nodes];
			for (int i = 0; i < g1->nr_nodes; ++i)
			{
				X[i] = new double[g1->nr_nodes];
				for (int j = 0; j < g1->nr_nodes; ++j)
					X[i][j] = (getSolution(varsx[i][j]) > 0.5) ? 1 : 0;
			}

			double **Y = new double*[g1->nr_nodes];
			for (int i = 0; i < g1->nr_nodes; ++i)
			{
				Y[i] = new double[g1->nr_nodes];
				for (int j = 0; j < g1->nr_nodes; ++j)
					Y[i][j] = (getSolution(varsy[i][j]) > 0.5) ? 1 : 0;
			}
			for (int i = 0; i < g1->nr_nodes; ++i)
			{
				int root;
				vector<bool> S(g1->nr_nodes, false);
				for (int q = 0; q < g1->nr_nodes; q++)
				{
					if (q == i)
						S[q] = true;
					if (X[i][q] == 1)
						root = q;
				}
				//do a BFS
				vector<int> children, parents;
				children.push_back(i);
				while(!children.empty())
				{ //do BFS
					parents = children;
					children.clear();
					for (int u = 0; u<parents.size(); ++u)
					{ //for each parent, examine the children
						int p = parents[u];
						//int j = -1;
						for (int v : g1->nb(p))
						{
							if (!S[v] && Y[v][p])
							{ //can only use vertices in S
								S[v] = true;
								children.push_back(v);
							}
						}
					}
				}
				if (S[root]) continue;
				GRBLinExpr expr1 = 0;
				GRBLinExpr expr2 = 0;
				for (int q = 0; q < g1->nr_nodes; q++)
				{
					if (S[q])
						expr1 += varsx[i][q];
				}
				for (int q = 0; q < g1->nr_nodes; q++)
				{
					if (!S[q]) continue;
					for (int u : g1->nb(q))
					{
						//int u = g1.adj[q][j];
						if (!S[u])
							expr2 += varsy[u][q];
					}
				}
				//Adding Lazycuts
				addLazy(expr1 + expr2 >= 1);
				numLazyCuts++;
			}

			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] X;
			delete[] Y;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}
