// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include "io.h"
#include "models.h"
#include <algorithm>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>

#ifndef sign
#define sign(x) (((x)>0)?1:((x)==0)?0:(-1))
#endif

#ifndef abs
#define abs(x) (((x)>0)?(x):(-(x)))
#endif

using namespace std;
//extern const char* gitversion;

#define DO_BATCH_OUTPUT

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(1);
    if (argc < 8) {
        printf("Usage: %s <dimacs> <distance> <population> <L|auto> <U|auto> <k> <model>\n\
  Available models:\n\
  \thess\t\tHess model\n\
  \tscf\t\tHess model with SCF\n\
  \tmcf1\t\tHess model with MCF1\n\
  \tmcf2\t\tHess model with MCF2\n\
  \tcut1\t\tHess model with CUT1\n\
  \tcut2\t\tHess modek with CUT2\n\
  \tul1\t\tU-L model type 1\n\
  \tul2\t\tU-L model type 2 (auto strong symmetry)\n", argv[0]);
        return 0;
    }
    // parse command line arguments
    char* dimacs_fname = argv[1];
    char* distance_fname = argv[2];
    char* population_fname = argv[3];
    int L = read_auto_int(argv[4], 0);
    int U = read_auto_int(argv[5], 0);

    // read inputs
    graph* g = 0;
    vector<vector<int> > dist;
    vector<int> population;
    if (read_input_data(dimacs_fname, distance_fname, population_fname, g, dist, population))
        return 1; // failure

    int k = read_auto_int("auto", g->get_k());

    if (L == 0 && U == 0)
        calculate_UL(population, k, &L, &U);

    //printf("Model input: L = %d, U = %d, k = %d\n", L, U, k);

    g->connect(dist);

    // check connectivity
    if (!g->is_connected())
    {
        printf("Problem is infeasible (not connected!)\n");
#ifdef DO_BATCH_OUTPUT
        printf("qwerky567: %s, disconnected\n", dimacs_fname);
#endif
        return 1;
    }

    if (g->nr_nodes <= 0)
    {
        fprintf(stderr, "empty graph\n");
        return 1;
    }

    if (dist.size() != g->nr_nodes || population.size() != g->nr_nodes)
    {
        fprintf(stderr, "dist/population size != n, expected %d\n", g->nr_nodes);
        return 1;
    }

    long total_pop = 0L;
    for (int p : population)
        total_pop += p;
    //printf("Model input: total population = %ld\n", total_pop);

    string arg_model = "hess";

    auto start = chrono::steady_clock::now();

    try {
        // initialize environment and create an empty model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // little ugly
        bool need_solution = true;
        GRBVar** x = 0;
        x = build_hess(&model, g, dist, population, L, U, k);

        auto start = chrono::steady_clock::now();
        //TODO change user-interactive?
        //model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour
        model.set(GRB_IntParam_Threads, 10);
        model.set(GRB_DoubleParam_NodefileStart, 10); // 10 GB
        model.set(GRB_IntParam_Method, 2);
        model.set(GRB_IntParam_OutputFlag, 0);
        //model.set(GRB_DoubleParam_MIPGap, 0);
        //model.set(GRB_DoubleParam_BarConvTol, 1e-12);
        //model.set(GRB_IntParam_CrossoverBasis, 1);
        model.set(GRB_IntParam_Crossover, 0);
        //optimize the model
        
        model.optimize();

        cout << fixed;

        GRBConstr *c = 0;

        c = model.getConstrs();

        int i;

        //print dual vars

        for (i = 0; i < g->nr_nodes; ++i)
        {
            cout << c[i].get(GRB_DoubleAttr_Pi) << endl;
        }

        vector<double> dualLambda(g->nr_nodes);

        for (; i < 2*g->nr_nodes; ++i)
        {
            cout << L*c[i].get(GRB_DoubleAttr_Pi) << endl;
        }


        for (; i < 3 * g->nr_nodes; ++i)
        {
            cout << U*c[i].get(GRB_DoubleAttr_Pi) << endl;
        }
        chrono::duration<double> duration = chrono::steady_clock::now() - start;
        cerr<<"Time elapsed: "<< duration.count() <<endl;
        //printf("Time elapsed: %lf seconds\n", duration.count()); // TODO use gurobi Runtime model attr
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (const char* msg) {
        cout << "Exception with message : " << msg << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }
    delete g;
    
    return 0;
}