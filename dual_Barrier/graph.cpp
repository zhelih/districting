#include "graph.h"
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <stack>
#include "rank.hpp"

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

    //printf("graph: %d nodes, %d edges\n", nr_nodes, nr_edges);

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

void graph::edgeClean(const vector<int>& population, int U)
{
    int numEdgeClean = 0;
    vector<int> nr_deleted(nr_nodes, 0);
    // O(m)
    for (int i = 0; i < nr_nodes; ++i)
        for (int j = 0; j < nb_[i].size(); ++j)
            if (population[i] + population[nb_[i][j]] > U)
            {
                numEdgeClean++;
                nb_[i][j] = -1;
                nr_deleted[i] += 1;
            }
    // now clean -1's from nb_
   // O(deleted)
    for (int i = 0; i < nr_nodes; ++i)
        if (nr_deleted[i])
        {
            int count = 0;
            for (int j = 0; j < nb_[i].size(); ++j)
                if (nb_[i][j] != -1)
                    swap(nb_[i][count++], nb_[i][j]);
            nb_[i].resize(nb_[i].size() - nr_deleted[i]);
        }
    cout << "# of cleaned edges: " << numEdgeClean / 2 << endl;
}

void graph::edgeCleanNeighbor(vector<int>& new_population, int s, int U)
{
    int numEdgeClean = 0;
    vector<int> nr_deleted(nr_nodes, 0);

    vector<bool> sClosedNeighbor(nr_nodes, false);
    sClosedNeighbor[s] = true;
    for (int i : nb(s))
        sClosedNeighbor[i] = true;

    for (int i = 0; i < nr_nodes; ++i)
    {
        if (sClosedNeighbor[i] == false) continue;
        for (int j = 0; j < nb_[i].size(); ++j)
        {
            if (sClosedNeighbor[nb_[i][j]] == false) continue;
            if (new_population[i] + new_population[nb_[i][j]] > U)
            {
                numEdgeClean++;
                nb_[i][j] = -1;
                nr_deleted[i] += 1;
            }
        }
    }
       
    for (int i = 0; i < nr_nodes; ++i)
        if (nr_deleted[i])
        {
            int count = 0;
            for (int j = 0; j < nb_[i].size(); ++j)
                if (nb_[i][j] != -1)
                    swap(nb_[i][count++], nb_[i][j]);
            nb_[i].resize(nb_[i].size() - nr_deleted[i]);
        }
    cout << "# of cleaned edges: " << numEdgeClean / 2 << endl;
}



vector<vector<int>> graph::FindBiconnectedComponents(vector<int> &AV, vector<bool> &deletedNodes)
{
    /* From Kgraph: I tried to use the naming conventions presented in Tarjan's 1972 paper */

    // declare vars
 
    vector< vector<int> > BC;		// biconnected components
    stack<int> le, re;				// used to store a stack of edges. le is "left edge" and re is "right edge". An edge is stored (le[i],re[i]). 

    //find connected components
    vector<vector<int>> connectedComponents = FindConnectedComponents(deletedNodes);

    //cerr << "# of connected components: " << connectedComponents.size() << endl;

    //for each connected component, find biconnected components
    for (int c = 0; c < connectedComponents.size(); ++c)
    {
        int v = connectedComponents[c][0];
        int u = -1, i = 0;
        vector<int> number(nr_nodes, (int)-1);
        vector<int> lowopt(nr_nodes, nr_nodes);

        // perform DFS-based algorithm
        Bico_Sub(v, u, i, number, lowopt, le, re, BC, deletedNodes);

        vector<int> countComp(nr_nodes, 0);
        for (int p = 0; p < BC.size(); ++p) // count how many components each vertex belongs to
            for (int q = 0; q < BC[p].size(); ++q)
                countComp[BC[p][q]]++;

        for (int p = 0; p < nr_nodes; ++p) // if a vertex belongs to >1 component, then it is a cut vertex
            AV[p] = countComp[p];
    }
    return BC;
}

void graph::Bico_Sub(int v, int u, int &i, vector<int> &number, vector<int> &lowopt, stack<int> &le, stack<int> &re, vector< vector<int>> &BC, vector<bool> &deletedNodes)
{
    i++;
    number[v] = i;
    lowopt[v] = number[v];
    int w;
    for (int w : nb(v))
    {
        if (deletedNodes[v] || deletedNodes[w]) continue;
        if (number[w] == -1)
        {
            le.push(v);
            re.push(w);
            Bico_Sub(w, v, i, number, lowopt, le, re, BC, deletedNodes);
            lowopt[v] = (int)min(lowopt[v], lowopt[w]);
            if (lowopt[w] >= number[v])
            {
                vector<int> temp_BC;
                vector<bool> bBC(nr_nodes, false);
                while (!le.empty() && !re.empty() && number[le.top()] >= number[w])
                {
                    bBC[le.top()] = true;
                    bBC[re.top()] = true;
                    le.pop();
                    re.pop();
                }
                if (!le.empty() && le.top() == v)
                {
                    bBC[le.top()] = true;
                    bBC[re.top()] = true;
                    le.pop();
                    re.pop();
                }
                else
                {
                    cout << "ERROR: edge (v,w) not on top of stack" << endl;
                }
                for (int p = 0; p < nr_nodes; ++p)
                    if (bBC[p])
                        temp_BC.push_back(p);
                BC.push_back(temp_BC);
            }
        }
        else if (number[w] < number[v] && w != u)
        {
            le.push(v);
            re.push(w);
            lowopt[v] = min(lowopt[v], number[w]);
        }
    }
}


void graph::connect(const vector<vector<int>>& dist)
{

    struct t_edge {
        int v1_graph;
        int v2_graph;

        int v1_comp;
        int v2_comp;

        int dist;
    };

    // run DFS to find connected components
    vector<int> comp(nr_nodes); // [i] component
    int nr_comp = 0; // number of connected components
    stack<int> s; // stack for the DFS
    vector<int> visited(nr_nodes, 0);
    for (int i = 0; i < nr_nodes; ++i) // start DFS
    {
        if (!visited[i])
        {
            s.push(i);
            while (!s.empty())
            {
                int v = s.top(); s.pop();
                if (!visited[v])
                {
                    visited[v] = true;
                    comp[v] = nr_comp;
                    for (int nb_v : nb(v))
                    {
                        if (!visited[nb_v])
                            s.push(nb_v);
                    }
                }
            }
            nr_comp++;
        }
    }

    if (nr_comp != 1)
    {
        //printf("Input graph is disconnected; adding these edges to make it connected:");
        vector<t_edge*> edges;
        // construct components reverse map for calculating distances
        vector<vector<int>> comp_rev(nr_comp);
        for (int i = 0; i < nr_nodes; ++i)
            comp_rev[comp[i]].push_back(i);
        // create the edge set

        // for every two components
        for (int comp1 = 0; comp1 < nr_comp; ++comp1)
            for (int comp2 = comp1 + 1; comp2 < nr_comp; ++comp2)
            {
                int c1_v_min = comp_rev[comp1][0]; int c2_v_min = comp_rev[comp2][0]; int d_min = dist[c1_v_min][c2_v_min];
                for (int c1_v : comp_rev[comp1])
                    for (int c2_v : comp_rev[comp2])
                        if (dist[c1_v][c2_v] < d_min)
                        {
                            c1_v_min = c1_v;
                            c2_v_min = c2_v;
                            d_min = dist[c1_v][c2_v];
                        }

                t_edge* edge = new t_edge;
                edge->v1_graph = c1_v_min;
                edge->v2_graph = c2_v_min;
                edge->v1_comp = comp1;
                edge->v2_comp = comp2;
                edge->dist = d_min;
                edges.push_back(edge);
            }
        // union-find for components
        union_of_sets u;
        u.dad = new int[nr_comp];
        u.rank = new int[nr_comp];
        for (int i = 0; i < nr_comp; ++i)
            MakeSet(u, i);

        // kruskal
        sort(edges.begin(), edges.end(), [](t_edge* e1, t_edge* e2) { return e1->dist < e2->dist; });
        for (t_edge* e : edges)
        {
            int r1 = Find(u, e->v1_comp);
            int r2 = Find(u, e->v2_comp);
            if (r1 != r2)
            {
                Union(u, r1, r2);
                add_edge(e->v1_graph, e->v2_graph);
            }
        }
        delete[] u.dad;
        delete[] u.rank;
        for (t_edge* e : edges)
            delete e;
    }
}

vector<vector<int>> graph::FindConnectedComponents(vector<bool> &deletedNodes)
{
    vector<vector<int>> connectedComponents;
    vector<bool> reached(nr_nodes, false);
    for (int i = 0; i < nr_nodes; ++i)
    {
        if (reached[i] || deletedNodes[i]) continue;
        //do a DFS to find nodes that are reachable from i
        reached[i] = true;
        vector<int> children;
        vector<int> oneComponent;
        oneComponent.push_back(i);
        children.push_back(i);
        while (!children.empty())
        { //do DFS
            int cur = children.back(); children.pop_back();
            for (int nb_cur : nb(cur))
            {
                if (!reached[nb_cur] && !deletedNodes[nb_cur])
                {
                    reached[nb_cur] = true;
                    children.push_back(nb_cur);
                    oneComponent.push_back(nb_cur);
                }
            }
        }
        connectedComponents.push_back(oneComponent);
    }
    return connectedComponents;
}

vector<int> graph::ShortestPathsNodeWeighted(int origin, const vector<int>& population, vector<bool> &deletedNodes, int total_pop)
{
    int v, u;
    int infty = total_pop + 1;
    vector<int> dist(nr_nodes, infty); //shortest distance from origin node to each other node. dist[i] = \infty means i not reachable
    vector<bool> Q(nr_nodes, true); // the set of vertices whose distance labels are not yet permanent
    dist[origin] = population[origin]; //the distance from the origin node is equal to its population
    int minDistance;
    int minVertex;
    //make Q[v]= false for deleted nodes
    //int numOfDeleted = 0;
    for (int i = 0; i < nr_nodes; ++i)
    {
        if (!deletedNodes[i]) continue;
        //dist[i] = infty;
        Q[i] = false;
        //numOfDeleted++;
    }
    
    // do node-weighted version of Dijkstra's shortest path algorithm.
    for (int i = 0; i < nr_nodes; ++i)
    {
        // find a vertex u from Q of minimum distance label dist[u]
        minDistance = infty;
        minVertex = -1;
        for (int v = 0; v < nr_nodes; ++v)
        {
            if (!Q[v]) continue;
            if (dist[v] < minDistance)
            {
                minDistance = dist[v];
                minVertex = v;
            }
        }
        //cerr << minVertex << endl;
        if (minVertex == -1) continue;
        // remove minVertex from Q
        Q[minVertex] = false;

        // update distance labels of neighbors
        for (int v : nb(minVertex))
        {
            if (Q[v] && dist[minVertex] + population[v] < dist[v])
            {
                dist[v] = dist[minVertex] + population[v];
            }
        }
    }
    return dist;
}