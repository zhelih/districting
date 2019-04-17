// source file for preprocess algorithm
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"

using namespace std;

vector<vector<int>> preprocess(graph* g, vector<int>& new_population, int L, int U, const vector<int>& population)
{
    vector<int> stem(g->nr_nodes);

    //clean edges of the input graph G (step 1)
    g->edgeClean(population, U);

    //initialization of stem (step 2)
    for (int i = 0; i < g->nr_nodes; ++i)
        stem[i] = i;

    //duplicate the input graph g (step 3)
    graph* g1 = g->duplicate();

    //compute the biconnected components (step 4) 
    vector<vector<int>> biconnectedComponents;
    vector<int> AV(g->nr_nodes);
    int cur;
    vector<bool> deletedNodes(g->nr_nodes, false);

    do {
        biconnectedComponents = g1->FindBiconnectedComponents(AV, deletedNodes);
        cur = FindMergableBiconnectedComponent(biconnectedComponents, new_population, population, AV, L);
        if (cur == -1) break;
        //update the population of stem (steps 7 and 10)
        int s;
        int curNode;
        int bcc_population = 0;
        for (int i = 0; i < biconnectedComponents[cur].size(); ++i)
        {
            curNode = biconnectedComponents[cur][i];
            if (AV[curNode] > 1)
                s = curNode;
            else
            {
                bcc_population += new_population[curNode];
                new_population[curNode] = 0;
                deletedNodes[curNode] = true;
            }
        }
        new_population[s] += bcc_population;

        //update the vector stem (steps 8 and 9)
        for (int i = 0; i < biconnectedComponents[cur].size(); ++i)
        {
            int node = biconnectedComponents[cur][i];
            stem[node] = s;
        }

        // cleaning in steps 11-13
        g->edgeCleanNeighbor(new_population, s, U);
        g1->edgeCleanNeighbor(new_population, s, U);

    } while (cur != -1);

    //cleaning steps 15 and 16
    QuickTestForInfeasibility(g1, new_population, deletedNodes, L, U);

    cerr << "# of merged nodes: " << count(deletedNodes.begin(), deletedNodes.end(), true) << endl;

    delete g1;

    return FindClustersFromStemVector(g, stem);
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
    //print clusters
    cerr << "Printing clusters obtained by function FindClustersFromStemVector" << endl;
    for (int i = 0; i < clusters.size(); i++)
    {
        cerr << "This is cluster: " << i << endl;
        for (int j = 0; j < clusters[i].size(); j++)
        {
            cerr << clusters[i][j] << endl;
        }
    }
    delete g2;
    return clusters;
}

void QuickTestForInfeasibility(graph* g, vector<int>& new_population, vector<bool>& deleted, int L, int U)
{
    //check overt feasibility
    for (int v = 0; v < g->nr_nodes; ++v)
    {
        if (new_population[v] > U)
        {
            printf("The instance is overt infeasible!\n");
            exit(0);
        }
    }

    //find underpopulated connected components in the auxiliary graph
    vector<bool> visited(g->nr_nodes, false);
    int it = 0;
    for (int i = 0; i < g->nr_nodes; ++i)
    {
        if (deleted[i] || visited[i]) continue;
        vector<int> children;
        int totalPopulation = new_population[i];
        children.push_back(i);
        while (!children.empty())
        { //do DFS
            int cur = children.back(); children.pop_back();
            for (int nb_cur : g->nb(cur))
            {
                if (!visited[nb_cur] && !deleted[nb_cur])
                {
                    visited[nb_cur] = true;
                    children.push_back(nb_cur);
                    totalPopulation += new_population[nb_cur];
                }
            }
        }
        if (totalPopulation < L)
        {
            printf("The instance is infeasible by an underpopulated component!\n");
            exit(0);
        }
    }
}

int FindMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, int L)
{
    int s = -1;
    for (int i = 0; i < biconnectedComponents.size(); ++i)
    {
        //if (!activeBiconnectedComponents[i]) continue;
        int cutVertexCounter = 0;
        int cutVertex;
        int bcc_population = 0;
        for (int j = 0; j < biconnectedComponents[i].size(); ++j)
        {
            int cur = biconnectedComponents[i][j];
            if (AV[cur] != 1)
            {
                ++cutVertexCounter;
                cutVertex = cur;
            }
            bcc_population += new_population[cur];
        }
        if (cutVertexCounter == 1 && bcc_population < L)
        {
            s = i;
            break;
        }
    }
    return s;
}
