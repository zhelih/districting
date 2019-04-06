// source file for lagrangian relaxation model
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "gurobi_c++.h"
#include "Models.h"

using namespace std;


vector<bool> solveInnerProblem(graph* g, double* x, vector<vector<bool>>& F_0, vector<vector<bool>>& F_1, int L, int U, int k, const vector<vector<int>>& clusters, const vector<vector<double>>& w, vector<vector<double>>& w_hat, vector<double>& W, const vector<int>& population)
{
    vector<bool> S(g->nr_nodes, false);
    // pass just "x" (the lagrangian concatenation of (alpha, lambda, upsilon) and inside of function create pointers to the subarrays, like Eugene mentioned.
    double *alpha = x;
    double *lambda = x + g->nr_nodes;
    double *upsilon = x + 2 * g->nr_nodes;
    // don't create w_hat from scratch each time. allocating a new set of memory for it each time will be time-consuming. Pass it as an argument?

    

    for (int i = 0; i < g->nr_nodes; i++)
    {
        double pOverL = static_cast<double>(population[i]) / static_cast<double>(L);
        double pOverU = static_cast<double>(population[i]) / static_cast<double>(U);
        for (int j = 0; j < g->nr_nodes; j++)
        {
            if (i == j)
                w_hat[i][j] = w[i][j] - alpha[i] - lambda[j] * pOverL + upsilon[j] * pOverU + lambda[j] - upsilon[j];
            else
            {
                w_hat[i][j] = w[i][j] - alpha[i] - lambda[j] * pOverL + upsilon[j] * pOverU;
            }
        }
    }


    //compute W_j for each vertex j \in V
    for (int c = 0; c < clusters.size(); c++)
    {
        for (int v = 0; v < clusters[c].size(); v++)
        {
            int j = clusters[c][v];
            //sum up \hat{w}_ij in the first term for all i \in \hat C 
            for (int u = 0; u < clusters[c].size(); u++)
            {
                int i = clusters[c][u];
                if (F_1[i][j]) continue;
                W[j] += w_hat[i][j];
            }
            //sum up \hat{w}_ij in the first term for all (i,j )\in F_1
            for (int i = 0; i < g->nr_nodes; i++)
            {
                if (F_1[i][j])
                    W[j] += w_hat[i][j];
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
                if (sum > 0)
                    W[j] += sum;
            }
        }
    }


    //find number of fixed centers 
    int fixed = 0;

    for (int i = 0; i < g->nr_nodes; i++)
    {
        for (int j = 0; j < g->nr_nodes; j++)
        {
            if (i != j) continue;
            if (F_1[i][i])
                fixed++;
        }
    }


    //solve the reduced Lagrangian
    vector<double> W_prime(g->nr_nodes, 0);
    W_prime = W;

    sort(W_prime.begin(), W_prime.end());

    vector<bool> selected(g->nr_nodes, false);
    int counter = 0;
    for (int i = 0; i < k - fixed; i++)
    {
        for (int j = 0; j < g->nr_nodes; j++)
        {
            if (W_prime[i] == W[j] && !selected[j] && counter < k - fixed)
            {
                selected[j] = true;
                S[j] = true;
                counter++;
            }
        }
    }
    return S;
}
