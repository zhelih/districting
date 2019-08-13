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

    //for (int i = 0; i < nr_nodes; ++i)
    //    cerr <<"vertex " <<i<<" : " << comp[i] << endl;

    fprintf(stderr, "nr_comp = %d\n", nr_comp);

    if (nr_comp == 1)
        printf("Graph is connected.\n");
    else
    {
        printf("Input graph is disconnected; adding these edges to make it connected:");
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
                printf(" { %d, %d }", e->v1_graph, e->v2_graph);
            }
        }
        printf("\n");

        delete[] u.dad;
        delete[] u.rank;
        for (t_edge* e : edges)
            delete e;
    }
}

int graph::get_k() const
{
  if(!(k > 0))
  {
    fprintf(stderr, "Value k is uninitialized!\n");
    exit(1);
  };
  return k;
}
