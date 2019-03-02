#include "graph.h"
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

//#define FROM_1

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

graph* from_dimacs(const char* fname)
{
    FILE* f = fopen(fname, "r");
    if (!f)
    {
        fprintf(stderr, "Cannot open file %s\n", fname);
        return 0;
    }

    char buf[1024]; int nr_nodes = 0, nr_edges = 0;

    // search for "n nodes edges"
    while (NULL != fgets(buf, sizeof(buf), f))
    {
        if (buf[0] == 'p')
        {
            sscanf(buf, "p edge %d", &nr_nodes);
            //      sscanf(buf, "p edge %d %d", &nr_nodes, &nr_edges);
            break;
        }
    }

    if (nr_nodes == 0 && nr_edges == 0)
    {
        fclose(f);
        fprintf(stderr, "Cannot found dimacs metadata in %s\n", fname);
        return 0;
    }

    //  fprintf(stderr, "Found metadata in %s : (%d, %d)\n", fname, nr_nodes, nr_edges);

    graph* g = new graph(nr_nodes);

    rewind(f);

    //read the edges
    while (NULL != fgets(buf, sizeof(buf), f))
    {
        if (buf[0] == 'e')
        {
            int from, to; double weight;
            sscanf(buf, "e %d %d %lf", &from, &to, &weight);
#ifdef FROM_1
            from--; to--;
#endif
            g->add_edge(from, to);//, weight);
        }
        else if (buf[0] == 'c' && buf[2] == 'k')
        {
            int k;
            sscanf(buf, "c k %d", &k);
            g->set_k(k);
        }
    }

    printf("graph: %d nodes, %d edges\n", nr_nodes, nr_edges);

    for (int i = 0; i < nr_nodes; ++i)
        sort(g->nb(i).begin(), g->nb(i).end());

    fclose(f);
    return g;
}

graph::graph(uint n) : k(0), nr_nodes(n)
{
    nb_.resize(n);
}

graph::~graph()
{
}

void graph::add_edge(uint i, uint j)
{
    if (is_edge(i, j))
        return;
    nb_[i].push_back(j);
    nb_[j].push_back(i);
}

void graph::remove_edge(uint i, uint j)
{
    if (!is_edge(i, j))
        return;

    int it = 0;
    for (int v : nb(i))
    {
        if (v == j)
            break;
        it++;
    }
    nb(i).erase(nb(i).begin() + it);

    it = 0;
    for (int v : nb(j))
    {
        if (v == i)
            break;
        it++;
    }
    nb(j).erase(nb(j).begin() + it);
}

bool graph::is_connected()
{
    vector<int> s;
    s.push_back(0);
    int *visited = new int[nr_nodes];
    for (uint i = 0; i < nr_nodes; ++i)
        visited[i] = 0;
    while (!s.empty())
    {
        int cur = s.back(); s.pop_back();
        visited[cur] = 1;
        auto nbs = nb(cur);
        for (auto it = nbs.begin(); it != nbs.end(); ++it)
        {
            if (!visited[*it])
            {
                s.push_back(*it);
            }
        }
    }
    bool res = true;
    for (unsigned int i = 0; i < nr_nodes; ++i)
        if (!visited[i])
        {
            res = false;
            break;
        }
    delete[] visited;
    return res;
}

bool graph::is_edge(uint i, uint j)
{
    for (uint k : nb_[i])
        if (k == j)
            return true;
    for (uint k : nb_[j])
        if (k == i)
            return true;
    return false;
}

void graph::edgeClean(vector<int>& population, int U)
{
    int numEdgeClean = 0;
    //remove unnecessary edges in input graph G
    for (int i = 0; i < nr_nodes; i++)
    {
        bool applyBreak = false;
        for (int j : nb(i))
        {
            if (population[i] + population[j] > U)
            {
                remove_edge(i, j);
                numEdgeClean++;
                applyBreak = true;
                break;
            }
        }
        if (applyBreak == true)
            i--;
    }
    cout << "# of cleaned edges: " << numEdgeClean << endl;
}

void graph::clean(vector<int>& new_population, vector<bool>& deleted, int L, int U, int& numOfEdgeDel, int& numOfNodeMerge)
{
    //check overt feasibility
    for (int v = 0; v < nr_nodes; v++)
    {
        if (new_population[v] > U)
        {
            printf ("The instance is overt infeasible!\n");
            exit(0);
        }
    }
    
    //remove unnecessary edges in the auxiliary graph
    for (int i = 0; i < nr_nodes; i++)
    {
        bool applyBreak = false;
        for (int j : nb(i))
        {
            if (new_population[i] + new_population[j] > U)
            {
                remove_edge(i, j);
                numOfEdgeDel++;
                applyBreak = true;
                break;
            }
        }
        if (applyBreak == true)
            i--;
    }

    //find underpopulated connected components in the auxiliary graph
    vector<bool> visited(nr_nodes, false);
    int it = 0;
    for (int i = 0; i < nr_nodes; i++)
    {
        if (deleted[i] == true)
            continue;
        if (visited[i] == true)
            continue;
        vector<int> children;
        int totalPopulation = new_population[i];
        children.push_back(i);
        while (!children.empty())
        { //do DFS
            int cur = children.back(); children.pop_back();
            for (int nb_cur : nb(cur))
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
            cout << "numOfNodeMerge: " << numOfNodeMerge << endl;
            cout << "numOfEdgeDel: " << numOfEdgeDel << endl;
            printf("The instance is infeasible by an underpopulated component!\n");
            exit(0);
        }
    }
}

vector<int> graph::findUnderPopulatedLeaves(vector<int> new_population, vector<bool> deleted, int L)
{
    vector<int> leaves;
    for (int i = 0; i <nr_nodes; i++)
    {
        if (deleted[i] == true)
            continue;
        int numNeigh = 0;
        for (int j : nb(i))
        {
            if (deleted[j] == true)
                continue;
            numNeigh++;
        }
        if (numNeigh == 1 && new_population[i] < L)
        {
            leaves.push_back(i);
        }
    }
    return leaves;
}

vector<int> graph::preprocess(vector<int>& new_population, vector<bool>& deleted, int L, int U)
{
    int numOfEdgeDel = 0;
    int numOfNodeMerge = 0;
    clean(new_population, deleted, L, U, numOfEdgeDel, numOfNodeMerge);
    vector<int> stem(nr_nodes,-1);
    vector<int> leaves = findUnderPopulatedLeaves (new_population, deleted, L);
    while (!leaves.empty())
    {
        int cur = leaves.back(); leaves.pop_back();
        int potStem;
        for (int j : nb(cur))
        {
           if (deleted[j] == true)
                continue;
           potStem = j;
        }
        stem[cur] = potStem;
        new_population[potStem] += new_population[cur];
        deleted[cur] = true;
        numOfNodeMerge++;
        remove_edge(cur, potStem);
        numOfEdgeDel++;
        clean(new_population, deleted, L, U, numOfEdgeDel, numOfNodeMerge);
        leaves = findUnderPopulatedLeaves(new_population, deleted, L);
    }
    cout << "numOfNodeMerge: " << numOfNodeMerge << endl;
    cout << "numOfEdgeDel: " << numOfEdgeDel << endl;
    return stem;
}
