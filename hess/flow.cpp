// source file for single and multi commodity flow formulations
#include <unordered_map>
#include <vector>
#include "gurobi_c++.h"
#include "graph.h"

void build_scf(GRBModel* model, GRBVar** x, graph* g, vector<int> stem)
{
    // create n^2 variables for arcs, presolve will eliminate unused
    int n = g->nr_nodes;

    // strengthening by stem inequalities
    for (int v = 0; v < n; ++v)
    {
        for (int i = 0; i < n; ++i)
        {
            if (stem[i] != -1)
                model->addConstr(x[i][v] - x[stem[i]][v] == 0);
        }
    }

    GRBVar** y = new GRBVar*[n];
    for (int i = 0; i < n; ++i)
        y[i] = model->addVars(n, GRB_BINARY);
    GRBVar** f = new GRBVar*[n];
    for (int i = 0; i < n; ++i)
        f[i] = model->addVars(n, GRB_CONTINUOUS);

    model->update();

    // add constraint (16b)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr expr = 0;
        for (int j : g->nb(i))
        {
            expr += y[j][i];
        }
        model->addConstr(expr + x[i][i] == 1);
    }

    // add constraint (16c)
    for (int i = 0; i < n; ++i)
    {
        GRBLinExpr expr = 0;
        for (int j : g->nb(i))
        {
            expr += f[i][j]; // in d_+
            expr -= f[j][i]; // in d_-
        }

        for (int j = 0; j < n; ++j)
            expr -= x[j][i];

        model->addConstr(expr == -1);
    }

    // add constraint (16d)
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
        {
            model->addConstr(f[i][j] - y[i][j] >= 0);
            model->addConstr(f[i][j] - n * y[i][j] <= 0);
        }

    // add constraint (Austin)
    for (int v = 0; v < n; ++v)
    {
        // for every edge here
        for (int i = 0; i < n; ++i)
            for (int j : g->nb(i))
                model->addConstr(x[i][v] + y[i][j] - x[j][v] <= 1);
    }
    model->update();

    model->write("debug_scf.lp");
}

void build_mcf0(GRBModel* model, GRBVar** x, graph* g, vector<int> stem)
{
    int n = g->nr_nodes;

    // strengthening by stem inequalities
    for (int v = 0; v < n; ++v)
    {
        for (int i = 0; i < n; ++i)
        {
            if (stem[i] != -1)
                model->addConstr(x[i][v] - x[stem[i]][v] == 0);
        }
    }

    // step 1 : hash edge (i,j) to i*n+j = h1
    // step 2 : hash every h1 to a number, resulting in exactly |E| variables

    std::unordered_map<int, int> hash_edges;
    int cur = 0;
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            hash_edges.insert(std::make_pair(i*n + j, cur++));

    // get number of edges
    int nr_edges = hash_edges.size();

    // add flow variables f[v][i,j]
    GRBVar**f = new GRBVar*[n]; // commodity type, v
    for (int v = 0; v < n; ++v)
        f[v] = model->addVars(nr_edges, GRB_CONTINUOUS); // the edge

      // add constraint (16b)
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            if (i == j) continue;
            GRBLinExpr expr = 0;
            for (int nb_i : g->nb(i))
            {
                expr += f[j][hash_edges[n*nb_i + i]]; // in d^- : edge (nb_i -- i)
                expr -= f[j][hash_edges[n*i + nb_i]]; // in d^+ : edge (i -- nb_i)
            }
            model->addConstr(expr == x[i][j]);
        }
    }

    // add constraint (16c)
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            if (i == j) continue;
            GRBLinExpr expr = 0;
            for (int nb_i : g->nb(i))
                expr += f[j][hash_edges[n*nb_i + i]]; // in d^- : edge (nb_i -- i)
            model->addConstr(expr <= (n - 1) * x[i][j]);
        }
    }

    // add constraint (16d) -- actually just fix individual vars to zero
    for (int j = 0; j < n; ++j)
        for (int i : g->nb(j))
            f[j][hash_edges[n*i + j]].set(GRB_DoubleAttr_UB, 0.); // in d^+ : edge (nb_j -- j)
}

void build_mcf1(GRBModel* model, GRBVar** x, graph* g, vector<int> stem)
{
    int n = g->nr_nodes;

    // strengthening by stem inequalities
    for (int v = 0; v < n; ++v)
    {
        for (int i = 0; i < n; ++i)
        {
            if (stem[i] != -1)
                model->addConstr(x[i][v] - x[stem[i]][v] == 0);
        }
    }

    // step 1 : hash edge (i,j) to i*n+j = h1
    // step 2 : hash every h1 to a number, resulting in exactly |E| variables

    std::unordered_map<int, int> hash_edges;
    int cur = 0;
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            hash_edges.insert(std::make_pair(i*n + j, cur++));

    // get number of edges
    int nr_edges = hash_edges.size();

    // add flow variables f[v][i,j]
    GRBVar**f = new GRBVar*[n]; // commodity type, v
    for (int v = 0; v < n; ++v)
        f[v] = model->addVars(nr_edges, GRB_CONTINUOUS); // the edge

                                 // add constraint (16b)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j) continue;
            GRBLinExpr expr = 0;
            for (int nb_j : g->nb(j))
            {
                expr += f[i][hash_edges[n*j + nb_j]]; // in d^+ : edge (j -- nb_j)
                expr -= f[i][hash_edges[n*nb_j + j]]; // in d^- : edge (nb_j -- j)
            }
            model->addConstr(expr == x[i][j]);
        }
    }

    // add constraint (16c) -- actually just fix individual vars to zero
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            f[i][hash_edges[n*i + j]].set(GRB_DoubleAttr_UB, 0.); // in d^+ : edge (i -- nb_i)

                                       // add constraint (16d)
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            for (int v = 0; v < n; ++v)
            {
                if (v == i || v == j) continue;
                model->addConstr(f[v][hash_edges[n*i + j]] <= f[j][hash_edges[n*i + j]]);
            }

    // add constraint (16e) 
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            f[j][hash_edges[n*i + j]].set(GRB_CharAttr_VType, GRB_BINARY);
}

void build_mcf2(GRBModel* model, GRBVar** x, graph* g, vector<int> stem)
{
    int n = g->nr_nodes;

    // strengthening by stem inequalities
    for (int v = 0; v < n; ++v)
    {
        for (int i = 0; i < n; ++i)
        {
            if (stem[i] != -1)
                model->addConstr(x[i][v] - x[stem[i]][v] == 0);
        }
    }

    // hash edges similarly to mcf1
    std::unordered_map<int, int> hash_edges;
    int cur = 0;
    for (int i = 0; i < n; ++i)
        for (int j : g->nb(i))
            hash_edges.insert(std::make_pair(i*n + j, cur++));

    // get number of edges
    int nr_edges = hash_edges.size();

    GRBVar ***f = new GRBVar**[n]; // f[ b ][ (i,j) ][ a ]
    for (int i = 0; i < n; ++i)
    {
        f[i] = new GRBVar*[nr_edges];
        for (int j = 0; j < nr_edges; ++j)
            f[i][j] = model->addVars(n - g->nb(i).size() - 1, GRB_CONTINUOUS); // V = { i } u N(i) u (V \ N[i])
    }

    // preprocess V \ N[i] sets
    std::vector<std::vector<int>> non_nbs(n);

    for (int i = 0; i < n; ++i)
    {
        std::vector<bool> nb(n, false);
        nb[i] = true;
        for (int j : g->nb(i))
            nb[j] = true;
        nb.flip();
        for (int j = 0; j < n; ++j)
            if (nb[j])
                non_nbs[i].push_back(j);

        // check
        if (non_nbs[i].size() != n - 1 - g->nb(i).size())
            throw "Internal Error : non nb size for mcf2";
    }

    // add constraint (19b)
    for (int b = 0; b < n; ++b)
        for (int a_i = 0; a_i < n - 1 - g->nb(b).size(); ++a_i)
        {
            GRBLinExpr expr = 0;
            for (int j : g->nb(b))
            {
                expr += f[b][hash_edges[b*n + j]][a_i]; // b -- j in d^+(b)
                expr -= f[b][hash_edges[j*n + b]][a_i]; // j -- b in d^-(b)
            }
            model->addConstr(expr == x[non_nbs[b][a_i]][b]);
        }

    // add constraint (19c)
    for (int b = 0; b < n; ++b)
        for (int a_i = 0; a_i < n - 1 - g->nb(b).size(); ++a_i)
            for (int i = 0; i < n; ++i)
            {
                if (i == non_nbs[b][a_i] || i == b) continue;
                GRBLinExpr expr = 0;
                for (int j : g->nb(i))
                {
                    expr += f[b][hash_edges[i*n + j]][a_i];
                    expr -= f[b][hash_edges[j*n + i]][a_i];
                }
                model->addConstr(expr == 0);
            }

    // add constraint (19d) -- actually just fix each var UB to zero
    for (int b = 0; b < n; ++b)
        for (int a_i = 0; a_i < n - 1 - g->nb(b).size(); ++a_i)
            for (int j : g->nb(b))
                f[b][hash_edges[j*n + b]][a_i].set(GRB_DoubleAttr_UB, 0.); // j -- b

                                              // add constraint (19e)
    for (int b = 0; b < n; ++b)
        for (int a_i = 0; a_i < n - 1 - g->nb(b).size(); ++a_i)
            for (int j = 0; j < n; ++j)
            {
                if (j == b) continue;
                GRBLinExpr expr = 0;
                for (int i : g->nb(j))
                    expr += f[b][hash_edges[i*n + j]][a_i]; // i -- j
                model->addConstr(expr <= x[j][b]);
            }
}
