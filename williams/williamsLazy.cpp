#include "LazyConstraints.h"
#include "gurobi_c++.h"
#include "KGraph.h"
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

long LazyConstraints::numCallbacks = 0;
double LazyConstraints::TotalCallbackTime = 0;
long LazyConstraints::numLazyCuts = 0;


GRBVar **vars;
GRBVar *vars2;
KGraph g1Cop;
KGraph d1;

long K1;
long L1;
long U1;
LazyConstraints::LazyConstraints(GRBVar **xvars, GRBVar *rvars, KGraph &g, KGraph &d, long &K, long &L, long &U)
{
	vars = xvars;
	vars2 = rvars;
	g1Cop.Duplicate(g);
	d1 = d;
	K1 = K;
	L1 = L;
	U1 = U;
};

void LazyConstraints::callback() {
	try {
		if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all vertex cut constraints so far.
		{
			//cerr << "hi" << endl;
			numCallbacks++;
			time_t start = clock();

			vector<bool> selectedEdges(g1Cop.m, false);
			vector<bool> selectedRoots(g1Cop.n, false);
			//vector<double> *x = new vector<double>[g1Cop.m];

			x = new long*[g1Cop.m];
			for (long i = 0; i < g1Cop.m; ++i)
				x[i] = new long[2];
			
			for (long i = 0; i < g1Cop.m; i++)
			{
				for (long j = 0; j < 2; j++)
				{
					if(getSolution(vars[i][j]) > 0.5)
						selectedEdges[i] = true;
				}
			}
			
			

			r = new long[g1Cop.n];
			for (long i = 0; i < g1Cop.n; i++)
			{
				if (getSolution(vars2[i]) > 0.5)
					selectedRoots[i] = true;
			}



			//cerr << "hi" << endl;

			//Run BFS
			long u,v,e,cutEdge;
			vector<long> labelVertex(g1Cop.n, -1);
			//vector<bool> markedEdge(g1Cop.m);
			vector<long> children;
			vector<long> parents;
			vector<long> cut;
			long sum =0;
			for (long i = 0; i < g1Cop.n; i++) //find root nodes
			{
				children.clear();
				parents.clear();
				cut.clear();
				sum = 0;
				if (selectedRoots[i])
				{
					//if (numCallbacks == 1)
					//	cerr << "Root: "<< i << endl;
					//cerr << d1.pop[i] << endl;
					sum += d1.pop[i];
					//cerr << i << endl;
					labelVertex[i] = i;
					//u = i;
					children.push_back(i);
					while (!children.empty())
					{
						parents = children;
						children.clear();
						for (long j = 0; j<parents.size(); j++)
						{
							u = parents[j];
							//cerr << "vertex " << u << ":" << endl;
							for (long aux = 0; aux < g1Cop.edj[u].size(); aux++)
							{
								e = g1Cop.edj[u][aux];
								//cerr << "edge: " << u << endl;
								if (selectedEdges[e]) 
								{
									for (long k = 0; k < 2; k++)
									{
										v = g1Cop.edge[e][k];
										if (labelVertex[v] != -1)
											continue;
										labelVertex[v] = i;
										//if (numCallbacks == 1)
										//	cerr << v << endl;
										children.push_back(v);
										sum += d1.pop[v];
									}
								}
							}
						}
					}
						
					for (long k = 0; k < g1Cop.m; k++)
					{
						if (selectedEdges[k])
							continue;
						if (labelVertex[g1Cop.edge[k][0]] == i && labelVertex[g1Cop.edge[k][1]] == i)
							continue;
						if (labelVertex[g1Cop.edge[k][0]] == i || labelVertex[g1Cop.edge[k][1]] == i)
						{
							cut.push_back(k);
						}
					}
						
					GRBLinExpr expr = 0;
					if (sum < L1)
					{
						//cerr << "Here is under populated!!!!!!!" << endl;
						for (long aux = 0; aux < cut.size(); aux++)
						{
							cutEdge = cut[aux];
							expr += vars[cutEdge][0] + vars[cutEdge][1];
						}
						addLazy(expr >= 1);
						numLazyCuts++;
						//break;
					}
					if (sum > U1)
					{
						//cerr << "Here is over populated!!!!!!!" << endl;
						for (long aux = 0; aux < g1Cop.n; aux++)
						{
							if (labelVertex[aux] == i)
							{
								expr += vars2[aux];
							}
						}
						for (long aux = 0; aux < cut.size(); aux++)
						{
							cutEdge = cut[aux];
							if (labelVertex[g1Cop.edge[cutEdge][0]] == i)
								expr += vars[cutEdge][0];
							if (labelVertex[g1Cop.edge[cutEdge][1]] == i)
								expr += vars[cutEdge][1];
						}
						addLazy(expr >= ceil(sum/U1));
						numLazyCuts++;
						//break;
					}
				}
			}

			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
		//cerr << "# of lazy cuts = " << numLazyCuts<<endl;
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}
