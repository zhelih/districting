// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"

void build_cut1(GRBModel* model, GRBVar** x, graph* g)
{
	// create n^2 variables, and set UB=0
	int n = g->nr_nodes;
	
	GRBVar** y = new GRBVar*[n];
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
}

void build_cut2(GRBModel* model, GRBVar** x, graph* g)
{
	throw "Unimplemented!";
}
