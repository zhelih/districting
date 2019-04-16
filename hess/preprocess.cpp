// source file for preprocess algorithm
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "gurobi_c++.h"
#include "Models.h"

using namespace std;

int FindOneMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, int L)
{
    int S = -1;
    for (int i = 0; i < biconnectedComponents.size(); ++i)
    {
        //if (!activeBiconnectedComponents[i]) continue;
        int noneOneCounter = 0;
        int noneOne;
        int sum = 0;
        for (int j = 0; j < biconnectedComponents[i].size(); ++j)
        {
            if (AV[biconnectedComponents[i][j]] != 1)
            {
                ++noneOneCounter;
                noneOne = biconnectedComponents[i][j];
            }
            sum += new_population[biconnectedComponents[i][j]];
        }
        if (noneOneCounter == 1 && sum < L)
        {
            S = i;
            break;
        }
    }
    return S;
}

vector<vector<int>> preprocess(graph* g, vector<int>& new_population, int L, int U, const vector<int>& population)
{
    vector<vector<int>> clusters;
    vector<int> stem(g->nr_nodes);
    int numOfEdgeDel = 0;
    int numOfNodeMerge = 0;
    //clean edges of the input graph G (step 1)
    g->edgeClean(population, U);

    //initialization of stem (step 2)
    for (int i = 0; i < g->nr_nodes; ++i)
        stem[i] = i;

    //duplicate the input graph g (step 3)
    graph* g1 = 0;
    g1 = g->duplicate();

    //compute the biconnected components (step 4) 
    vector<vector<int>> biconnectedComponents;
    vector<int> AV(g->nr_nodes);
    int cur;
    vector<bool> deletedNodes(g->nr_nodes, false);
    biconnectedComponents = FindBiconnectedComponents(g1, AV, deletedNodes);

    do {
        cur = FindOneMergableBiconnectedComponent(biconnectedComponents, new_population, population, AV, L);
        if (cur >= 0)
        {
            //update the population of stem (steps 7 and 10)
            int s;
            int curNode;
            int sum = 0;
            for (int i = 0; i < biconnectedComponents[cur].size(); ++i)
            {
                curNode = biconnectedComponents[cur][i];
                if (AV[curNode] > 1)
                    s = curNode;
                else
                {
                    sum += new_population[curNode];
                    new_population[curNode] = 0;
                    deletedNodes[curNode] = true;
                }
            }
            new_population[s] += sum;

            //update the vector stem (steps 8 and 9)
            for (int i = 0; i < biconnectedComponents[cur].size(); ++i)
            {
                int node = biconnectedComponents[cur][i];
                if (node == s)
                    continue;
                stem[node] = s;
            }

            // cleaning in steps 11-13
            for (int v : g1->nb(s))
            {
                if (new_population[v] + new_population[s] > U)
                {
                    g1->remove_edge(v, s);
                    g->remove_edge(v, s);
                }
            }

            biconnectedComponents = FindBiconnectedComponents(g1, AV, deletedNodes);
            cur = FindOneMergableBiconnectedComponent(biconnectedComponents, new_population, population, AV, L);
        }
    } while (cur != -1);

    //translate stem to set family \mathcal{C} 
    clusters = FindClustersFromStemVector(g,stem);

    for (int i = 0; i < clusters.size(); i++)
    {
        cerr << "This is cluster: " << i << endl;
        for (int j = 0; j < clusters[i].size(); j++)
        {
            cerr << clusters[i][j] << endl;
        }
    }
    //cleaning steps 15 and 16
    g1->clean(new_population, deletedNodes, L, U, numOfEdgeDel, numOfNodeMerge);

    delete g1;

    return clusters;
}

vector<vector<int>> FindClustersFromStemVector(graph* g, vector<int>& stem)
{
    //translate stem to set family \mathcal{C} 
    vector<vector<int>> clusters;
    graph* g2 = new graph(g->nr_nodes);

    for (int i = 0; i < g2->nr_nodes; ++i)
        if (stem[i] != i)
            g2->nb(stem[i]).push_back(i);

    vector<bool> R(g2->nr_nodes, false);
    for (int i = 0; i < g2->nr_nodes; ++i)
    {
        if (stem[i] != i || R[i]) continue;
        //do a DFS to find nodes that are reachable from i
        R[i] = true;
        vector<int> children;
        vector<int> oneCluster;
        oneCluster.push_back(i);
        children.push_back(i);
        while (!children.empty())
        { //do DFS
            int cur = children.back(); children.pop_back();
            for (int nb_cur : g2->nb(cur))
            {
                if (!R[nb_cur])
                {
                    R[nb_cur] = true;
                    children.push_back(nb_cur);
                    oneCluster.push_back(nb_cur);
                }
            }
        }
        clusters.push_back(oneCluster);
    }
    delete g2;
    return clusters;
}
