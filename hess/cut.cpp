// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include <chrono>

class CutCallback : public HessCallback
{
    // memory for a callback
private:
    int* visited; // dfs marks
    int* aci; // A(C_i) set
    std::vector<int> s; // stack for DFS
    std::vector<int> cc; // connected component for a vertex
public:
    CutCallback(hess_params& p, graph *g_) : HessCallback(p, g_)
    {
        visited = new int[n];
        aci = new int[n];
        s.reserve(n);
        cc.reserve(n);
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
                        //run dfs from j and mark cc
                        cur_cc++; // we don't really need info about cc # but it is neat to have it
                        s.clear(); s.push_back(j); cc[j] = cur_cc;
                        cc_heads.push_back(j);
                        while(!s.empty())
                        {
                          int cur = s.back(); s.pop_back();
                          for(int nb_cur : g->nb(cur))
                            if(x_val[nb_cur][i] > 0.5 && !visited[nb_cur] && cc[j] == 0)
                            {
                              cc[nb_cur] = cur_cc;
                              s.push_back(nb_cur);
                            }
                        }
                      }
                    for(int j : cc_heads) 
                    {
                        // compute i-j separator, A(C_i) is already computed)
                        GRBLinExpr expr = 0;
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
                                    if (aci[nb_cur])
                                        expr += X(nb_cur, i);
                                    else s.push_back(nb_cur);
                                }
                            }
                        }
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

HessCallback* build_cut(GRBModel* model, hess_params& p, graph* g)
{
    model->getEnv().set(GRB_IntParam_LazyConstraints, 1); // turns off presolve!!!
    CutCallback* cb = new CutCallback(p, g);
    model->setCallback(cb);
    model->update();
    return cb;
}
