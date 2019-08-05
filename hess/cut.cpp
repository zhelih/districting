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
  int* aci; // A(C_b) set
  std::vector<int> s; // stack for DFS
  std::vector<int> cc; // connected component for a vertex
  std::vector<int> dist;
  bool is_lcut;
  int U;
public:
  CutCallback(hess_params& p, graph *g_, const vector<int>& pop_, bool is_lcut_, int U_) : HessCallback(p, g_, pop_), is_lcut(is_lcut_), U(U_)
  {
    visited = new int[n];
    aci = new int[n];
    s.reserve(n);
    cc.resize(n);
    dist.resize(n);
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
      for (int b = 0; b < n; ++b)
      {
        if (x_val[b][b] > 0.5) // b is a clusterhead
        {
          // run DFS from b on C_b, compute A(C_b) to save time later
          // clean A(C_b) and visited
          for (int j = 0; j < n; ++j)
          {
            aci[j] = 0;
            visited[j] = 0;
          }

          // run DFS on C_b
          s.clear(); s.push_back(b); visited[b] = 1;
          while (!s.empty())
          {
            int cur = s.back(); s.pop_back();
            for (int nb_cur : g->nb(cur))
              if (x_val[nb_cur][b] > 0.5) // if nb_cur is in C_b
              {
                if (!visited[nb_cur])
                {
                  visited[nb_cur] = 1;
                  s.push_back(nb_cur);
                }
              }
              else aci[nb_cur] = 1; // nb_cur is a neighbor of a vertex in C_b, thus in A(C_b)
          }

          // here if C_b is connected, all vertices in C_b must be visited
          // since we want to add cut for every connected component reamining there, we will mark cc's
          int cur_cc = 0;
          fill(cc.begin(), cc.end(), 0);
          vector<int> cc_heads;
          for (int j = 0; j < n; ++j)
            if (x_val[j][b] > 0.5 && !visited[j] && cc[j] == 0)
            {
              int cc_max_pop_node = j;
              //run dfs from j and mark cc
              cur_cc++; // we don't really need info about cc # but it is neat to have it
              s.clear(); s.push_back(j); cc[j] = cur_cc;
              while (!s.empty())
              {
                int cur = s.back(); s.pop_back();
                for (int nb_cur : g->nb(cur))
                  if (x_val[nb_cur][b] > 0.5 && !visited[nb_cur] && cc[nb_cur] == 0)
                  {
                    cc[nb_cur] = cur_cc;
                    s.push_back(nb_cur);
                    if (population[nb_cur] > population[cc_max_pop_node])
                      cc_max_pop_node = nb_cur;
                  }
              }
              cc_heads.push_back(cc_max_pop_node);
            }
          for (int a : cc_heads)
          {
            // compute i-j separator, A(C_b) is already computed)
            GRBLinExpr expr = 0;
            unordered_set<int> C;
            // start DFS from j to find R_i
            for (int k = 0; k < n; ++k)
              visited[k] = 0;
            s.clear(); s.push_back(a); visited[a] = 1;
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
            if (is_lcut)
            {
              // refine set C
              for (auto it_c = C.begin(); it_c != C.end(); )
              {
                int c = *it_c;
                // find distance from a to b through c but not remaining C
                // priority queue Dijkstra
                // priority queue pair is <weight, vertex>, min weight on top
                priority_queue< pair<int, int>, vector <pair<int, int>>, greater<pair<int, int>> > pq;
                fill(dist.begin(), dist.end(), INT_MAX);
                const vector<int>& p = population; // alias
                pq.push(make_pair(p[a], a));
                dist[a] = p[a];
                while (!pq.empty())
                {
                  int u = pq.top().second; pq.pop();
                  for (int nb_u : g->nb(u))
                  {
                    if ((nb_u == c || C.count(nb_u) == 0) && dist[nb_u] > dist[u] + p[nb_u])
                    {
                      dist[nb_u] = dist[u] + p[nb_u];
                      pq.push(make_pair(dist[nb_u], nb_u));
                    }
                  }
                }
                if (dist[b] > U)
                  it_c = C.erase(it_c);
                else
                  ++it_c;
              }
            }
            for (int c : C)
              expr += X(c, b);
            expr -= X(a, b); // RHS
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
  CutCallback* cb = new CutCallback(p, g, population, is_lcut, U);
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
  // fix x_ab=0 if dist_{G,p}(a,b)>U 
  vector<int> dist(g->nr_nodes);
  int countFixed = 0;
  for (int b = 0; b < g->nr_nodes; ++b)
  {
    // compute population-weighted distance from a to all other nodes.
    priority_queue< pair<int, int>, vector <pair<int, int>>, greater<pair<int, int>> > pq;
    fill(dist.begin(), dist.end(), INT_MAX);
    pq.push(make_pair(population[b], b));
    dist[b] = population[b];
    while (!pq.empty())
    {
      int u = pq.top().second; pq.pop();
      for (int nb_u : g->nb(u))
      {
        if (dist[nb_u] > dist[u] + population[nb_u])
        {
          dist[nb_u] = dist[u] + population[nb_u];
          pq.push(make_pair(dist[nb_u], nb_u));
        }
      }
    }
    for (int a = 0; a < g->nr_nodes; ++a)
      if (dist[a] > U && IS_X(a, b))
      {
        model->addConstr(X_V(a, b) == 0);
        countFixed++;
      }
  }
  cout << "Number of vars fixed in lcut initialization = " << countFixed << endl;

  return build_cut_(model, p, g, population, true, U);
}
