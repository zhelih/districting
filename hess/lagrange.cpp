// source file for lagrangian relaxation model
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"

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
    sort(W_indices.begin(), W_indices.end(), [&](int i1, int i2) { return W[i1] < W[i2]; });

    // null grad
    int it = 0;
    f_val = 0.;
    for(int it = 0; it < g->nr_nodes; ++it)
    {
      grad[it] = 1.;
      f_val += alpha[it];
    }
    for(; it < 3*g->nr_nodes; ++it)
      grad[it] = 0.;

    for (int i = 0; i < g->nr_nodes; ++i)
    {
        for (int j = 0; j < g->nr_nodes; ++j)
        {
            if (F_1[i][j])
            {
                // setting x[i][j] = true
                grad[i] -= 1;
                grad[j+g->nr_nodes] -= static_cast<double>(population[i]) / static_cast<double>(L);
                grad[j+2*g->nr_nodes] += static_cast<double>(population[i]) / static_cast<double>(U);
                f_val += w_hat[i][j];
            }
        }
    }

    for (int i = 0; i < k - fixed; i++)
    {
        // setting x[W_indices[i]][W_indices[i]] = true
        S[W_indices[i]] = true;
        grad[W_indices[i]+g->nr_nodes] += 1;
        grad[W_indices[i]+2*g->nr_nodes] -= 1;
        f_val += w_hat[W_indices[i]][W_indices[i]];
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
                if(F_1[i][j]) continue;
                // setting x[i][j] = true
                grad[i] -= 1;
                grad[j+g->nr_nodes] -= static_cast<double>(population[i]) / static_cast<double>(L);
                grad[j+2*g->nr_nodes] += static_cast<double>(population[i]) / static_cast<double>(U);
                f_val += w_hat[i][j];
            }
        }
    }
}
