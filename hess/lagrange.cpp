// source file for lagrangian relaxation model
#include <cstdio>
#include <cassert>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
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



void solveInnerProblem(graph* g, const double* multipliers, const vector<vector<bool>>& F_0, const vector<vector<bool>>& F_1,
    int L, int U, int k, const vector<vector<int>>& clusters, const vector<int>& population,
    const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, vector<bool>& S, double* grad, double& f_val)
{
    const double *alpha = multipliers;
    const double *lambda = multipliers + g->nr_nodes;
    const double *upsilon = multipliers + 2 * g->nr_nodes;

    for (int i = 0; i < g->nr_nodes; i++)
    {
        double pOverL = static_cast<double>(population[i]) / static_cast<double>(L);
        double pOverU = static_cast<double>(population[i]) / static_cast<double>(U);
        for (int j = 0; j < g->nr_nodes; j++)
        {
            if (i == j)
                w_hat[i][j] = w[i][j] - alpha[i] - abs(lambda[j]) * pOverL + abs(upsilon[j]) * pOverU + abs(lambda[j]) - abs(upsilon[j]);
            else
                w_hat[i][j] = w[i][j] - alpha[i] - abs(lambda[j]) * pOverL + abs(upsilon[j]) * pOverU;
        }
    }

    //compute W_j for each vertex j \in V
    for (int c = 0; c < clusters.size(); c++)
    {
        for (int v = 0; v < clusters[c].size(); v++)
        {
            int j = clusters[c][v];

            W[j] = 0;
            //sum up \hat{w}_ij in the first term
            for (int c1 = 0; c1 < clusters.size(); c1++)
            {
                if (c1 == c || F_1[clusters[c1][0]][j])
                {
                    for (int h = 0; h < clusters[c1].size(); h++)
                    {
                        int a = clusters[c1][h];
                        W[j] += w_hat[a][j];
                    }
                }
            }
            //sum up \hat{w}_ij in the second term
            for (int l = 0; l < clusters.size(); l++)
            {
                if (l == c) continue;
                double sum = 0;
                for (int h = 0; h < clusters[l].size(); h++)
                {
                    int a = clusters[l][h];
                    if (!F_0[a][j] && !F_1[a][j])
                        sum += w_hat[a][j];
                }
                W[j] += min(sum, 0);
            }
        }
    }

    //find number of fixed centers 
    int fixed = 0;

    for (int i = 0; i < g->nr_nodes; i++)
    {
        if (F_1[i][i])
            fixed++;
    }

    //solve the reduced Lagrangian
    vector<int> W_indices(W.size());
    for (size_t i = 0; i < W.size(); ++i)
        W_indices[i] = i;
    sort(W_indices.begin(), W_indices.end(), [&W](int i1, int i2) { return W[i1] < W[i2]; });

    // null grad
    int it = 0;
    f_val = 0.;
    for (int it = 0; it < g->nr_nodes; ++it)
    {
        grad[it] = 1.;
        f_val += alpha[it];
    }
    for (; it < 3 * g->nr_nodes; ++it)
        grad[it] = 0.;

    unordered_set<pair<int, int>, pair_hash> h;

    for (int i = 0; i < g->nr_nodes; ++i)
    {
        for (int j = 0; j < g->nr_nodes; ++j)
        {
            if (F_1[i][j])
            {
                if (i != j)
                {
                    // setting x[i][j] = true
                    grad[i] -= 1;
                    grad[j + g->nr_nodes] -= static_cast<double>(population[i]) / static_cast<double>(L);
                    grad[j + 2 * g->nr_nodes] += static_cast<double>(population[i]) / static_cast<double>(U);
                    h.insert(make_pair(i, j));
                }
                else
                {
                    h.insert(make_pair(i, i));
                    grad[i + g->nr_nodes] += 1;
                    grad[i + 2 * g->nr_nodes] -= 1;
                }
                f_val += w_hat[i][j];
            }
        }
    }

    for (int i = 0; i < k - fixed; i++)
    {
        // setting x[W_indices[i]][W_indices[i]] = true
        S[W_indices[i]] = true;
        grad[W_indices[i] + g->nr_nodes] += 1;
        grad[W_indices[i] + 2 * g->nr_nodes] -= 1;
        f_val += w_hat[W_indices[i]][W_indices[i]];
        h.insert(make_pair(W_indices[i], W_indices[i]));
    }

    // fix S[i][j] for all i \in C (i != j) if w_hat_{Cj} is nonnegative
    for (int j = 0; j < g->nr_nodes; ++j)
    {
        if (!S[j]) continue;
        for (int c = 0; c < clusters.size(); ++c)
        {
            double sum = 0;
            int i = -1;
            for (int b = 0; b < clusters[c].size(); ++b)
            {
                i = clusters[c][b];
                sum += w_hat[i][j];
            }
            if (sum > 0 || F_0[clusters[c][0]][j]) continue;
            for (int b = 0; b < clusters[c].size(); ++b)
            {
                i = clusters[c][b];
                if (F_1[i][j]) continue;
                // setting x[i][j] = true
                if (h.count(make_pair(i, j)) > 0)
                    continue;
                f_val += w_hat[i][j];
                if (i != j)
                {
                    grad[i] -= 1;
                    grad[j + g->nr_nodes] -= static_cast<double>(population[i]) / static_cast<double>(L);
                    grad[j + 2 * g->nr_nodes] += static_cast<double>(population[i]) / static_cast<double>(U);
                }
                else
                {
                    grad[i + g->nr_nodes] += 1;
                    grad[i + 2 * g->nr_nodes] -= 1;
                }
            }
        }
    }
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
        if (!S[j] && f_val + W[j] - maxW > UB) // is it possible to have x_jj=1 in an optimal solution?
            for (int i = 0; i < n; ++i) // the answer is no. Fix x_ij=0 for all i
                F_0[i][j] = true;
        else
        {
            for (int i = 0; i < n; ++i)
            {
                if (i == j) continue;
                if (!S[j] && f_val + W[j] - maxW + max(0, w_hat[i][j]) > UB) // is it possible to have x_ij=1 in an optimal solution?
                    F_0[i][j] = true;
                if (S[j] && f_val + max(0, w_hat[i][j]) > UB)
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
              if(F_0[i][j]) continue; // must skip from F_0
              if(F_1[i][j]) // must add from F_1
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
        if(F_0[i1][i1] + 2*F_1[i1][i1] == F_0[i2][i2] + 2*F_1[i2][i2]) // check if both are in the same set or neither
            return W[i1] < W[i2];
        if(F_0[i1][i1] || F_1[i2][i2]) return false;
        if(F_1[i1][i1] || F_0[i2][i2]) return true;
      return W[i1] < W[i2];// never reached
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
