// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include <chrono>
#include <unordered_set>
#include <queue> // priority queue

#ifndef INT_MAX
#include <limits>
#define INT_MAX (std::numeric_limits<int>::max())
#endif

class CutCallback : public HessCallback
{
    // memory for a callback
private:
    int* visited; // dfs marks
    int* aci; // A(C_i) set
    std::vector<int> s; // stack for DFS
    std::vector<int> cc; // connected component for a vertex
    bool is_lcut;
    int U;
public:
    CutCallback(hess_params& p, graph *g_, const vector<int>& pop_, bool is_lcut_, int U_) : HessCallback(p, g_, pop_), is_lcut(is_lcut_), U(U_)
    {
        visited = new int[n];
        aci = new int[n];
        s.reserve(n);
        cc.resize(n);
    }
    virtual ~CutCallback()
    {
        delete[] aci;
        delete[] visited;
    }
protected:
    void callback();
};

void CutCallback::callback()
{
    using namespace std;
    try
    {
        if (where == GRB_CB_MIPSOL)
        {
            ++numCallbacks;
            auto start = chrono::steady_clock::now();

            populate_x(); // from HessCallback

            // try clusterheads
            for (int i = 0; i < n; ++i)
            {
                if (x_val[i][i] > 0.5) // i is a clusterhead
                {
                    // run DFS from i on C_i, compute A(C_i) to save time later
                    // clean A(C_i) and visited
                    for (int j = 0; j < n; ++j)
                    {
                        aci[j] = 0;
                        visited[j] = 0;
                    }

                    // run DFS on C_i
                    s.clear(); s.push_back(i);
                    while (!s.empty())
                    {
                        int cur = s.back(); s.pop_back();
                        visited[cur] = 1;
                        for (int nb_cur : g->nb(cur))
                            if (x_val[nb_cur][i] > 0.5) // if nb_cur is in C_i
                            {
                                if (!visited[nb_cur])
                                    s.push_back(nb_cur);
                            }
                            else aci[nb_cur] = 1; // nb_cur is a neighbor of a vertex in C_i, thus in A(C_i)a
                    }

                    // here if C_i is connected, all vertices in C_i must be visited
                    // since we want to add cut for every connected component reamining there, we will mark cc's
                    int cur_cc = 0;
                    fill(cc.begin(), cc.end(), 0);
                    vector<int> cc_heads;
                    for(int j = 0; j < n; ++j)
                      if(x_val[j][i] > 0.5 && !visited[j] && cc[j] == 0)
                      {
                        int cc_max_pop_node = j;
                        //run dfs from j and mark cc
                        cur_cc++; // we don't really need info about cc # but it is neat to have it
                        s.clear(); s.push_back(j); cc[j] = cur_cc;
                        while(!s.empty())
                        {
                          int cur = s.back(); s.pop_back();
                          for(int nb_cur : g->nb(cur))
                            if(x_val[nb_cur][i] > 0.5 && !visited[nb_cur] && cc[nb_cur] == 0)
                            {
                              cc[nb_cur] = cur_cc;
                              s.push_back(nb_cur);
                              if(population[nb_cur] > population[cc_max_pop_node])
                                cc_max_pop_node = nb_cur;
                            }
                        }
                        cc_heads.push_back(cc_max_pop_node);
                      }
                    for(int j : cc_heads) 
                    {
                        // compute i-j separator, A(C_i) is already computed)
                        GRBLinExpr expr = 0;
                        unordered_set<int> C;
                        // start DFS from j to find R_i
                        for (int k = 0; k < n; ++k)
                            visited[k] = 0;
                        s.clear(); s.push_back(j); visited[j] = 1;
                        while (!s.empty())
                        {
                            int cur = s.back(); s.pop_back();
                            for (int nb_cur : g->nb(cur))
                            {
                                if (!visited[nb_cur])
                                {
                                    visited[nb_cur] = 1;
                                    if (aci[nb_cur]) C.insert(nb_cur); else s.push_back(nb_cur);
                                }
                            }
                        }
                        if(is_lcut)
                        {
                          // refine set C
                          for(int c : C)
                          {
                            // find distance from a to b through c but not remaining C
                            // priority queue Dijkstra
                            // priority queue pair is <weight, vertex>, min weight on top
                            priority_queue< pair<int,int>, vector <pair<int,int>> , greater<pair<int,int>> > pq;
                            vector<int> dist(n, INT_MAX);
                            const vector<int>& p = population; // alias
                            pq.push(make_pair(p[j], j));
                            dist[j] = p[j];
                            while(!pq.empty())
                            {
                              int u = pq.top().second; pq.pop();
                              for(int nb_u : g->nb(u))
                              {
                                if((nb_u == c || C.count(nb_u) == 0) && dist[nb_u] > dist[u] + p[nb_u])
                                {
                                  dist[nb_u] = dist[u] + p[nb_u];
                                  pq.push(make_pair(dist[nb_u], nb_u));
                                }
                              }
                            }
                            if(dist[i] > U)
                              C.erase(i);
                          }
                        }
                        for(int c : C)
                          expr += X(c, i);
                        expr -= X(j,i); // RHS
                        addLazy(expr >= 0);
                        ++numLazyCuts;
                    }
                }
            }
            chrono::duration<double> d = chrono::steady_clock::now() - start;
            callbackTime += d.count();
        }
    }
    catch (GRBException e)
    {
        fprintf(stderr, "Error number: %d\n", e.getErrorCode());
        fprintf(stderr, "%s\n", e.getMessage().c_str());
    }
    catch (...)
    {
        fprintf(stderr, "Error during callback\n");
    }
}

HessCallback* build_cut_(GRBModel* model, hess_params& p, graph* g, const vector<int>& population, bool is_lcut, int U)
{
    model->getEnv().set(GRB_IntParam_LazyConstraints, 1); // turns off presolve!!!
    CutCallback* cb = new CutCallback(p, g, population, is_lcut, 0);
    model->setCallback(cb);
    model->update();
    return cb;
}

HessCallback* build_cut(GRBModel* model, hess_params& p, graph* g, const vector<int>& population)
{
  return build_cut_(model, p, g, population, false, 0);
}
HessCallback* build_lcut(GRBModel* model, hess_params& p, graph* g, const vector<int>& population, int U)
{
  return build_cut_(model, p, g, population, true, U);
}
