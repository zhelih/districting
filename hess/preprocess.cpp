// source file for lagrangian relaxation model
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
    vector<vector<int>> clusters;
    vector<int> stem (g->nr_nodes);
    int numOfEdgeDel = 0;
    int numOfNodeMerge = 0;
    //clean edges of the input graph G (step 1)
    g->edgeClean(population, U);

    //initialization of stem (step 2)
    for (int i = 0; i < g->nr_nodes; i++)
        stem[i] = i;

    //duplicate the input graph g (step 3)
    graph* g1 = 0;
    g1 = g->duplicate();

    //compute the biconnected components (step 4) 
    vector<vector<int>> biconnectedComponents;
    vector<int> AV(g->nr_nodes);
    vector<int> mergableBiconnectedComponent;
    vector<bool> deletedNodes(g->nr_nodes, false);
    biconnectedComponents = FindBiconnectedComponents(g1, AV, deletedNodes);
    vector<bool> activeBiconnectedComponents(biconnectedComponents.size(), true);

    //find mergable biconnected components
    mergableBiconnectedComponent = FindMergableBiconnectedComponent(biconnectedComponents, new_population, population, AV, activeBiconnectedComponents, L);

    //step 5
    bool flag = !mergableBiconnectedComponent.empty();
    while (flag)
    {
        //select a mergable connected component
        int cur = mergableBiconnectedComponent.back();

        //update the population of stem (steps 7 and 10)
        int s;
        int curNode;
        int sum = 0;
        for (int i = 0; i < biconnectedComponents[cur].size(); i++)
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
        int counter = 0;
        for (int i = 0; i < biconnectedComponents[cur].size(); i++)
        {

            if (biconnectedComponents[cur][i] == s)
                continue;
            counter++;
            stem[biconnectedComponents[cur][i]] = s;
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

        for (int i = 0; i < biconnectedComponents.size(); i++)
            activeBiconnectedComponents[i] = true;

        for (int i = 0; i < biconnectedComponents.size(); i++)
        {
            int counter = 0;
            for (int j = 0; j < biconnectedComponents[i].size(); j++)
            {
                if (deletedNodes[biconnectedComponents[i][j]])
                    counter++;
            }
            if (counter == biconnectedComponents[i].size() - 1)
                activeBiconnectedComponents[i] = false;
        }
        mergableBiconnectedComponent = FindMergableBiconnectedComponent(biconnectedComponents, new_population, population, AV, activeBiconnectedComponents, L);
        flag = !mergableBiconnectedComponent.empty();
    }

    //translate stem to set family \mathcal{C} 
    graph* g2 = new graph(g->nr_nodes);

    for (int i = 0; i < g2->nr_nodes; i++)
        if (stem[i] != i)
            g2->add_edge(stem[i], i);

    vector<bool> R(g2->nr_nodes, false);
    for (int i = 0; i < g2->nr_nodes; i++)
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
    //cerr << "Size of cluster: " << clusters.size() << endl;

    for (int i = 0; i < g2->nr_nodes; i++)
    {
        cerr << "stem of " << i << " is " << stem[i] << endl;
    }

    cerr << "Size of Cluster: " << clusters.size() << endl;;

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

    return clusters;
}

vector<int> FindMergableBiconnectedComponent(vector<vector<int>>& biconnectedComponents, vector<int>& new_population, const vector<int>& population, vector<int>& AV, vector<bool>& activeBiconnectedComponents, int L)
{
    vector<int> S;
    for (int i = 0; i < biconnectedComponents.size(); i++)
    {
        if (!activeBiconnectedComponents[i]) continue;
        int noneOneCounter = 0;
        int noneOne;
        int sum = 0;
        for (int j = 0; j < biconnectedComponents[i].size(); j++)
        {
            if (AV[biconnectedComponents[i][j]] != 1)
            {
                noneOneCounter++;
                noneOne = biconnectedComponents[i][j];
            }
            sum += new_population[biconnectedComponents[i][j]];
        }
        if (noneOneCounter == 1 && sum < L)
            S.push_back(i);
    }
    return S;
}
