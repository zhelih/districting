// source file for lagrangian relaxation model
#include <cstdio>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "models.h"
#include "ralg/ralg.h"

struct pair_hash {
	inline std::size_t operator()(const std::pair<int, int> & v) const {
		return v.first * 31 + v.second;
	}
};
using namespace std;

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef abs
#define abs(x) (((x)>0)?(x):(-(x)))
#endif

double solveLagrangian(graph* g, const vector<vector<double>>& w, const vector<int> &population, int L, int U, int k, 
	vector<vector<double>>& LB0, vector<vector<double>>& LB1, vector<int> &lagrangianCenters)
{
	double LB = -INFINITY;

	vector<double> W(g->nr_nodes, 0);
	vector<vector<double>> w_hat(g->nr_nodes, vector<double>(g->nr_nodes));

	vector<bool> currentCenters(g->nr_nodes);	// centers from most recent inner problem
	vector<bool> bestCenters(g->nr_nodes);		// centers from inner problem with best LB

	int dim = 3 * g->nr_nodes;
	double * bestMultipliers = new double[dim]; 
	double * multipliers = new double[dim];

	auto cb_grad_func = [g, &w, &population, L, U, k, &W, &w_hat, &currentCenters, &bestCenters, &LB, &LB0, &LB1, dim](const double* multipliers, double& f_val, double* grad) 
	{
		solveInnerProblem(g, multipliers, L, U, k, population, w, w_hat, W, grad, f_val, currentCenters);
		update_LB0_and_LB1(W, currentCenters, f_val, w_hat, LB0, LB1);

		// update incubments?
		if (f_val > LB)
		{
			LB = f_val;
			bestCenters = currentCenters;
		}
		return true;
	};

	//auto start = chrono::steady_clock::now();

	// run ralg
	
	cerr << "\n\n**Currently incapable of using/writing hot start**\n\n";

	// try to load hot start if any
	//if (ralg_hot_start)
	//	read_ralg_hot_start(ralg_hot_start_fname, multipliers, dim);
	//else
	for (int i = 0; i < dim; ++i)
		multipliers[i] = 1.; // whatever
	
	ralg_options opt = defaultOptions; opt.output_iter = 1; opt.is_monotone = false;
	LB = ralg(&opt, cb_grad_func, dim, multipliers, bestMultipliers, RALG_MAX); // lower bound from lagrangian

	// dump result to "state_model.hot"
	/*string fn1 = string(dimacs_fname);
	string hsfn = fn1.substr(0, 2) + "_" + arg_model + ".hot";
	int len1 = hsfn.length() + 1;
	char* char_array1 = new char[len1];
	strcpy(char_array1, hsfn.c_str());
	dump_ralg_hot_start(char_array1, res, dim);*/

	int p = 0;
	for (int i = 0; i < g->nr_nodes; ++i)
		if (bestCenters[i])
		{
			lagrangianCenters[p] = i;
			p++;
		}

	delete[] multipliers;
	delete[] bestMultipliers;

	return LB;
}

void update_LB0_and_LB1(const vector<double>& W, const vector<bool>& currentCenters, double f_val, 
	const vector<vector<double>> &w_hat, vector< vector<double> > &LB0, vector< vector<double> > &LB1)
{
	int n = currentCenters.size();
	double maxW = -INFINITY;
	double minW = INFINITY;

	// determine values for minW and maxW
	for (int i = 0; i < n; ++i)
		if (currentCenters[i])
			maxW = max(maxW, W[i]);

	for (int i = 0; i < n; ++i)
		if (!currentCenters[i])
			minW = min(minW, W[i]);

	// update LB1
	for (int j = 0; j < n; ++j)
	{
		if (!currentCenters[j])
		{
			//if (maxW == -INFINITY) continue;
			LB1[j][j] = max(LB1[j][j], f_val + W[j] - maxW);
			for (int i = 0; i < n; ++i)
			{
				if (i == j) continue;
				LB1[i][j] = max(LB1[i][j], f_val + W[j] - maxW + max(0, w_hat[i][j]));
			}
		}
		else
		{
			for (int i = 0; i < n; ++i)
			{
				if (i == j) continue;
				LB1[i][j] = max(LB1[i][j], f_val + max(0, w_hat[i][j]));
			}
		}
	}

	// TODO currently not updating LB0, since x[i][j]=1 safe fixings are not common????
}

void lagrangianBasedSafeFixing(vector<vector<bool>>& F_0, vector<vector<bool>>& F_1,
	const vector<vector<int>>& clusters, vector<double>& W, const vector<bool>& S, const double f_val, const double UB, const vector<vector<double>> &w_hat)
{
	int n = S.size();
	double maxW = -INFINITY;
	double minW = INFINITY;

	// determine values for minW and maxW
	for (int i = 0; i < n; ++i)
		if (S[i] && !F_1[i][i])
			maxW = max(maxW, W[i]);

	for (int i = 0; i < n; ++i)
		if (!S[i] && !F_0[i][i])
			minW = min(minW, W[i]);

	//// fix some vars x_jj to zero?
	//for (int j = 0; j < n; ++j)
	//	if (!S[j] && f_val + W[j] - maxW > UB) // is it possible to have x_jj=1 in an optimal solution?
	//		for (int i = 0; i < n; ++i) // the answer is no. Fix x_ij=0 for all i
	//			F_0[i][j] = true;


	// fix some vars x_ij to zero?
	for (int j = 0; j < n; ++j)
	{
		if (!S[j])
		{
			if (f_val + W[j] - maxW > UB) // is it possible to have x_jj=1 in an optimal solution?)
			{
				for (int i = 0; i < n; ++i) // the answer is no. Fix x_ij=0 for all i
					F_0[i][j] = true;
			}
			else
			{
				for (int i = 0; i < n; ++i)
					if (i != j && f_val + W[j] - maxW + max(0, w_hat[i][j]) > UB)
						F_0[i][j] = true;
			}
		}
		else // S[j]
		{
			for (int i = 0; i < n; ++i)
			{
				if (i != j && f_val + max(0, w_hat[i][j]) > UB)
					F_0[i][j] = true;
			}
		}
	}

	// fix some vars x_jj to one?
	for (int c = 0; c < clusters.size(); ++c)
	{
		for (int p = 0; p < clusters[c].size(); ++p)
		{
			int j = clusters[c][p];

			// is it possible to have x_jj=0 in an optimal solution?
			if (S[j] && f_val - W[j] + minW > UB)
			{
				// the answer is no. Fix x_ij=1 for all i in the same cluster as j
				for (int p1 = 0; p1 < clusters[c].size(); ++p1)
				{
					int i = clusters[c][p1];
					F_1[i][j] = true;
					for (int j1 = 0; j1 < n; ++j1)
						if (j1 != j)
							F_0[i][j1] = true;
				}
			}
		}
	}

	// report the number of fixings
	int numFixedZero = 0;
	int numFixedOne = 0;
	int numUnfixed = 0;
	int numCentersLeft = 0;
	for (int i = 0; i < n; ++i)
	{
		if (!F_0[i][i]) numCentersLeft++;
		for (int j = 0; j < n; ++j)
		{
			if (F_0[i][j]) numFixedZero++;
			else if (F_1[i][j]) numFixedOne++;
			else numUnfixed++;
		}
	}
	cerr << endl;
	cerr << "Number of variables fixed to zero = " << numFixedZero << endl;
	cerr << "Number of variables fixed to one  = " << numFixedOne << endl;
	cerr << "Number of variables not fixed     = " << numUnfixed << endl;
	cerr << "Number of centers left            = " << numCentersLeft << endl;
	cerr << "Percentage of vars fixed = " << (double)(n*n - numUnfixed) / (n*n) << endl;
}

void solveInnerProblem(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
	const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& currentCenters)
{
	const double *alpha = multipliers;
	const double *lambda = multipliers + g->nr_nodes;
	const double *upsilon = multipliers + 2 * g->nr_nodes;

	// reset centers
	for (int i = 0; i < g->nr_nodes; ++i)
		currentCenters[i] = false;

	// update w_hat
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		double pOverL = static_cast<double>(population[i]) / static_cast<double>(L);
		double pOverU = static_cast<double>(population[i]) / static_cast<double>(U);
		for (int j = 0; j < g->nr_nodes; ++j)
		{
			// w[i][i] == 0?
			w_hat[i][j] = w[i][j] - alpha[i] - abs(lambda[j]) * pOverL + abs(upsilon[j]) * pOverU;
			if (i == j)
				w_hat[i][j] += abs(lambda[j]) - abs(upsilon[j]);
		}
	}

	// recompute W_j, the minimum obj value for district centered at j
	for (int j = 0; j < g->nr_nodes; ++j)
	{
		W[j] = w_hat[j][j]; 
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i != j && w_hat[i][j] < 0) 
				W[j] += w_hat[i][j];
	}

	// select k smallest
	vector<int> W_indices(W.size());
	for (size_t i = 0; i < W.size(); ++i)
		W_indices[i] = i;
	std::sort(W_indices.begin(), W_indices.end(), [&W](int i1, int i2) {
		return W[i1] < W[i2];
	});

	// compute f_val
	f_val = 0.;
	for (int i = 0; i < g->nr_nodes; ++i)
		f_val += alpha[i];

	for (int j = 0; j < k; ++j)
	{
		int v = W_indices[j];
		f_val += W[v];
		currentCenters[v] = true;
	}

	// compute grad
	// A
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		grad[i] = 1.;
		for (int j = 0; j < k; ++j)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[i] -= 1;
	}
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		grad[i + g->nr_nodes] = 0.;
		grad[i + 2 * g->nr_nodes] = 0.;
	}
	// L
	for (int j = 0; j < k; ++j)
	{
		grad[g->nr_nodes + W_indices[j]] = 1;
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[g->nr_nodes + W_indices[j]] -= static_cast<double>(population[i]) / static_cast<double>(L);
	}
	//U
	for (int j = 0; j < k; ++j)
	{
		grad[2 * g->nr_nodes + W_indices[j]] = -1;
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[2 * g->nr_nodes + W_indices[j]] += static_cast<double>(population[i]) / static_cast<double>(U);
	}

	// signify the gradient
	// L
	for (int i = 0; i < g->nr_nodes; ++i)
		if (lambda[i] < 0)
			grad[i + g->nr_nodes] = -grad[i + g->nr_nodes];
	// U
	for (int i = 0; i < g->nr_nodes; ++i)
		if (upsilon[i] < 0)
			grad[i + 2 * g->nr_nodes] = -grad[i + 2 * g->nr_nodes];
}

void eugene_inner(graph* g, const double* multipliers, int L, int U, int k, const vector<int>& population,
	const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, double* grad, double& f_val, vector<bool>& S,
	const vector<vector<bool>>& F_0, const vector<vector<bool>>& F_1)
{
	const double *alpha = multipliers;
	const double *lambda = multipliers + g->nr_nodes;
	const double *upsilon = multipliers + 2 * g->nr_nodes;

	// TODO compute S only when needed
	for (int i = 0; i < g->nr_nodes; ++i)
		S[i] = false;

	for (int i = 0; i < g->nr_nodes; ++i)
	{
		double pOverL = static_cast<double>(population[i]) / static_cast<double>(L);
		double pOverU = static_cast<double>(population[i]) / static_cast<double>(U);
		for (int j = 0; j < g->nr_nodes; ++j)
		{
			// w[i][i] == 0?
			w_hat[i][j] = w[i][j] - alpha[i] - abs(lambda[j]) * pOverL + abs(upsilon[j]) * pOverU;
			if (i == j)
				w_hat[i][j] += abs(lambda[j]) - abs(upsilon[j]);
		}
	}

	// W_j describes the minimum value what happpens if j is a clusterhead
	for (int j = 0; j < g->nr_nodes; ++j)
	{
		W[j] = w_hat[j][j]; // weight of clusterhead node
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i != j)
			{
				if (F_0[i][j]) continue; // must skip from F_0
				if (F_1[i][j]) // must add from F_1
				{
					W[j] += w_hat[i][j];
					continue;
				}
				// else
				W[j] += min(0., w_hat[i][j]); // add all negative weights from nodes
			}
	}

	// select k smallest
	// must include for F_1[i][i]
	// must skip for F_0[i][i]
	// assuming at most k F_1[i][i] are set to true
	vector<int> W_indices(W.size());
	for (size_t i = 0; i < W.size(); ++i)
		W_indices[i] = i;
	sort(W_indices.begin(), W_indices.end(), [&W, &F_0, &F_1](int i1, int i2) {
		if (F_0[i1][i1] + 2 * F_1[i1][i1] == F_0[i2][i2] + 2 * F_1[i2][i2]) // check if both are in the same set or neither
			return W[i1] < W[i2];
		return F_1[i1][i1] || F_0[i2][i2];
	});

	// compute f_val
	f_val = 0.;
	for (int i = 0; i < g->nr_nodes; ++i)
		f_val += alpha[i];
	for (int j = 0; j < k; ++j)
	{
		f_val += W[W_indices[j]];
		S[W_indices[j]] = true;
	}

	// compute grad
	// A
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		grad[i] = 1.;
		for (int j = 0; j < k; ++j)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[i] -= 1;
	}
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		grad[i + g->nr_nodes] = 0.;
		grad[i + 2 * g->nr_nodes] = 0.;
	}
	// L
	for (int j = 0; j < k; ++j)
	{
		grad[g->nr_nodes + W_indices[j]] = 1;
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[g->nr_nodes + W_indices[j]] -= static_cast<double>(population[i]) / static_cast<double>(L);
	}
	//U
	for (int j = 0; j < k; ++j)
	{
		grad[2 * g->nr_nodes + W_indices[j]] = -1;
		for (int i = 0; i < g->nr_nodes; ++i)
			if (i == W_indices[j] || w_hat[i][W_indices[j]] < 0)
				grad[2 * g->nr_nodes + W_indices[j]] += static_cast<double>(population[i]) / static_cast<double>(U);
	}

	// signify the gradient
	// L
	for (int i = 0; i < g->nr_nodes; ++i)
		if (lambda[i] < 0)
			grad[i + g->nr_nodes] = -grad[i + g->nr_nodes];
	// U
	for (int i = 0; i < g->nr_nodes; ++i)
		if (upsilon[i] < 0)
			grad[i + 2 * g->nr_nodes] = -grad[i + 2 * g->nr_nodes];
}
