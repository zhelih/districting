// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"

void build_cut1(GRBModel* model, GRBVar** x, graph* g)
{
	// create variables for arcs
	int n = g->nr_nodes;

	GRBVar** y = new GRBVar*[n];
	for (int i = 0; i < n; ++i)
		y[i] = model->addVars(n, GRB_BINARY);

	// add constraint (23b)
	for (int i = 0; i < n; ++i)
	{
		GRBLinExpr expr = 0;
		for (int j : g->nb(i))
		{
			expr += y[j][i];
		}
		model->addConstr(expr + x[i][i] == 1);
	}
}

void build_cut2(GRBModel* model, GRBVar** x, graph* g)
{
	throw "Unimplemented!";
}
