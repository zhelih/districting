// source file for cut based formulations
#include "graph.h"
#include "gurobi_c++.h"
#include "models.h"
#include <chrono>

//Lazy constraints for solving CUT1 in political redistricting
class Cut1Callback : public HessCallback
{
private:
    GRBVar **grb_y;
    double **y;
public:
    Cut1Callback(GRBVar **grb_x, GRBVar **yvars, graph *g, cvv& F0, cvv& F1) : HessCallback(grb_x, g, F0, F1), grb_y(yvars)
    {
        y = new double*[n];
        for (int i = 0; i < n; ++i)
            y[i] = new double[n];
    }
    virtual ~Cut1Callback()
    {
        for (int i = 0; i < n; ++i)
            delete[] y[i];
        delete[] y;
    }
protected:
    void callback();
    void populate_y()
    {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                y[i][j] = getSolution(grb_y[i][j]);
    }
};

void Cut1Callback::callback()
{
    using namespace std;
    try
    {
        if (where == GRB_CB_MIPSOL) // Found an integer ``solution'' that satisfies all cut constraints so far.
        {
            ++numCallbacks;
            auto start = chrono::steady_clock::now();

            populate_x(); // from HessCallback 
            populate_y(); // from Cut1Callback

            for (int i = 0; i < n; ++i)
            {
                int root;
                vector<bool> R(n, false);
                R[i] = true;
                for (int j = 0; j < n; ++j)
                {
                    if (x_val[i][j] > 0.5)
                        root = j;
                }
                //do a DFS to find nodes that can reach i
                vector<int> children;
                children.push_back(i);
                while (!children.empty())
                { //do DFS
                    int cur = children.back(); children.pop_back();
                    for (int nb_cur : g->nb(cur))
                    {
                        if (!R[nb_cur] && y[nb_cur][cur] + y[cur][nb_cur] > 0.5)
                        { //can only use vertices in S
                            R[nb_cur] = true;
                            children.push_back(nb_cur);
                        }
                    }
                }
                if (R[root]) continue;
                GRBLinExpr expr = 0;
                for (int j = 0; j < n; ++j)
                {
                    if (R[j])
                    {
                        expr += X(i,j);
                        for (int u : g->nb(j))
                        {
                            if (!R[u])
                                expr += grb_y[u][j];
                        }
                    }
                }
                //Adding Lazycuts
                addLazy(expr >= 1);
                ++numLazyCuts;
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

HessCallback* build_cut1(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1)
{
    // create n^2 variables, and set UB=0
    model->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    int n = g->nr_nodes;

    GRBVar** y = new GRBVar*[n]; //FIXME ever deleted?
    for (int i = 0; i < n; ++i)
    {
        GRBVar *y_temp = new GRBVar[n];
        for (int j = 0; j < n; ++j)
        {
            // set UBs to zero for all y vars
            y_temp[j] = model->addVar(0.0, 0.0, 0.0, GRB_BINARY);
        }
        y[i] = y_temp;
    }

    // add constraint (23b)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr expr = 0;
        for (int j : g->nb(i))
        {
            // set UBs to one for all arcs
            y[j][i].set(GRB_DoubleAttr_UB, 1.0);
            expr += y[j][i];
        }
        model->addConstr(expr + X(i,i) == 1);
    }
    //FIXME delete
    Cut1Callback* cb = new Cut1Callback(x, y, g, F0, F1); // tell Gurobi which function generates the lazy cuts.
    model->setCallback(cb);
    model->update();

    return cb;
}

class Cut2Callback : public HessCallback
{
    // memory for a callback
private:
    int* visited; // dfs marks
    int* aci; // A(C_i) set
    std::vector<int> s; // stack for DFS
public:
    Cut2Callback(GRBVar** grb_x_, graph *g_, cvv& F0, cvv& F1) : HessCallback(grb_x_, g_, F0, F1)
    {
        visited = new int[n];
        aci = new int[n];
        s.reserve(n);
    }
    virtual ~Cut2Callback()
    {
        delete[] aci;
        delete[] visited;
    }
protected:
    void callback();
};

void Cut2Callback::callback()
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
            bool done = false;
            for (int i = 0; i < n && !done; ++i)
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
                    for (int j = 0; j < n; ++j)
                        if (x_val[j][i] > 0.5 && !visited[j])
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
                            done = true;
                            break;
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

HessCallback* build_cut2(GRBModel* model, GRBVar** x, graph* g, cvv& F0, cvv& F1)
{
    model->getEnv().set(GRB_IntParam_LazyConstraints, 1); // turns off presolve!!!
    Cut2Callback* cb = new Cut2Callback(x, g, F0, F1);
    model->setCallback(cb);
    model->update();
    return cb;
}
