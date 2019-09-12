// source file for lagrangian relaxation model
#include <cstdio>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <queue>
#include "common.h"
#include "graph.h"
#include "models.h"
#include "ralg/ralg.h"
#include "io.h"

double solveLagrangian(graph* g, const vector<vector<double>>& w, const vector<int> &population, int L, int U, int k, 
  vector<vector<double>>& LB1, bool ralg_hot_start, const char* ralg_hot_start_fname, const run_params& rp, bool exploit_contiguity)
{
  double LB = -MYINFINITY;

  vector<double> W(g->nr_nodes, 0);
  vector<vector<double>> w_hat(g->nr_nodes, vector<double>(g->nr_nodes));

  vector<bool> currentCenters(g->nr_nodes); // centers from most recent inner problem

  int dim = 3 * g->nr_nodes;
  double * bestMultipliers = new double[dim]; 
  double * multipliers = new double[dim];

  auto cb_grad_func = [g, &w, &population, L, U, k, &W, &w_hat, &currentCenters, &LB, &LB1, dim, exploit_contiguity](const double* multipliers, double& f_val, double* grad) 
  {
    solveInnerProblem(g, multipliers, L, U, k, population, w, w_hat, W, grad, f_val, currentCenters);
    if (exploit_contiguity)
      update_LB_contiguity(g, W, currentCenters, f_val, w_hat, LB1);
    else
      update_LB(W, currentCenters, f_val, w_hat, LB1);

    // update incubments?
    if (f_val > LB)
      LB = f_val;

    return true;
  };

  // try to load hot start if any
  if (ralg_hot_start)
    read_ralg_hot_start(ralg_hot_start_fname, multipliers, dim);
  else
    for (int i = 0; i < dim; ++i)
      multipliers[i] = 1.; // whatever

  ralg_options opt = defaultOptions; opt.output_iter = 1; opt.is_monotone = false;
  if (ralg_hot_start) opt.itermax = 100;
  LB = ralg(&opt, cb_grad_func, dim, multipliers, bestMultipliers, RALG_MAX); // lower bound from lagrangian

  // dump result to "state_model.hot"
  dump_ralg_hot_start(rp, bestMultipliers, dim, LB);

  delete [] multipliers;
  delete [] bestMultipliers;

  return LB;
}

void update_LB(const vector<double>& W, const vector<bool>& currentCenters, double f_val, 
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB1)
{
  int n = currentCenters.size();
  double maxW = -MYINFINITY;
  double minW = MYINFINITY;

  // determine values for minW and maxW
  for (int i = 0; i < n; ++i)
    if (currentCenters[i])
      maxW = mymax(maxW, W[i]);

  for (int i = 0; i < n; ++i)
    if (!currentCenters[i])
      minW = mymin(minW, W[i]);

  // update LB1
  for (int j = 0; j < n; ++j)
  {
    if (!currentCenters[j])
    {
      //if (maxW == -MYINFINITY) continue;
      LB1[j][j] = mymax(LB1[j][j], f_val + W[j] - maxW);
      for (int i = 0; i < n; ++i)
      {
        if (i == j) continue;
        LB1[i][j] = mymax(LB1[i][j], f_val + W[j] - maxW + mymax(0, w_hat[i][j]));
      }
    }
    else
    {
      for (int i = 0; i < n; ++i)
      {
        if (i == j) continue;
        LB1[i][j] = mymax(LB1[i][j], f_val + mymax(0, w_hat[i][j]));
      }
    }
  }
}

void update_LB_contiguity(graph* g, const vector<double>& W, const vector<bool>& currentCenters, double f_val,
  const vector<vector<double>> &w_hat, vector< vector<double> > &LB1)
{
  int n = currentCenters.size();
  double maxW = -MYINFINITY;

  // determine value for maxW
  for (int i = 0; i < n; ++i)
    if (currentCenters[i])
      maxW = mymax(maxW, W[i]);

  // compute special distances 
  vector<double> dist(g->nr_nodes);
  for (int j = 0; j < n; ++j)
  {
    // a particular shortest path computation from j to all nodes
    priority_queue< pair<double, int>, vector <pair<double, int>>, greater<pair<double, int>> > pq;
    for (int i = 0; i < n; ++i) dist[i] = DBL_MAX;
    pq.push(make_pair(0., j)); // copy constructor?
    dist[j] = 0.; // NB: not zero here!

    while (!pq.empty())
    {
      int u = pq.top().second;
      pq.pop();
      for (int nb : g->nb(u)) {
        int weight = mymax(0, w_hat[nb][j]);
        if (dist[nb] > dist[u] + weight) {
          dist[nb] = dist[u] + weight;
          pq.push(make_pair(dist[nb], nb));
        }
      }
    }

    // update LB1[][]
    if (currentCenters[j])
    {
      for (int i = 0; i < n; ++i)
        LB1[i][j] = mymax(LB1[i][j], f_val + dist[i]);
    }
    else
    {
      for (int i = 0; i < n; ++i)
        LB1[i][j] = mymax(LB1[i][j], f_val - maxW + W[j] + dist[i]);
    }
  }
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
      w_hat[i][j] = w[i][j] - alpha[i] - myabs(lambda[j]) * pOverL + myabs(upsilon[j]) * pOverU;
      if (i == j)
        w_hat[i][j] += myabs(lambda[j]) - myabs(upsilon[j]);
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
