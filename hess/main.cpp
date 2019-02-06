// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include "io.h"
#include "models.h"
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>

using namespace std;
extern const char* gitversion;

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(1);
    printf("Districting, build %s\n", gitversion);
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

    int k = read_auto_int(argv[6], g->get_k());

    if (L == 0 && U == 0)
        calculate_UL(population, k, &L, &U);

    printf("Model input: L = %d, U = %d, k = %d\n", L, U, k);

    // check connectivity
    if (!g->is_connected())
    {
        printf("Problem is infeasible (not connected!)\n");
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
    printf("Model input: total population = %ld\n", total_pop);

    string arg_model = argv[argc - 1];


    try {
        // initialize environment and create an empty model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // little ugly
        bool need_solution = true;
        GRBVar** x = 0;
        if (arg_model != "ul1" && arg_model != "ul2")
            x = build_hess(&model, g, dist, population, L, U, k);

        HessCallback* cb = 0;

        if (arg_model == "scf")
            build_scf(&model, x, g);
        else if (arg_model == "mcf1")
            build_mcf1(&model, x, g);
        else if (arg_model == "mcf2")
            build_mcf2(&model, x, g);
        else if (arg_model == "cut1")
            cb = build_cut1(&model, x, g);
        else if (arg_model == "cut2")
            cb = build_cut2(&model, x, g);
        else if (arg_model == "shirabe")
            build_shirabe(&model, x, g);
        else if (arg_model == "ul1") {
            x = build_UL_1(&model, g, population, k);
            need_solution = false;
        }
        else if (arg_model == "ul2") {
            x = build_UL_2(&model, g, population, k);
            need_solution = false;
        }
        else if (arg_model != "hess") {
            fprintf(stderr, "ERROR: Unknown model : %s\n", arg_model.c_str());
            exit(1);
        }

        //TODO change user-interactive?
        model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour
        model.set(GRB_IntParam_Threads, 10);
        model.set(GRB_DoubleParam_NodefileStart, 10); // 10 GB

                                                      //optimize the model
        auto start = chrono::steady_clock::now();
        model.optimize();
        chrono::duration<double> duration = chrono::steady_clock::now() - start;
        printf("Time elapsed: %lf seconds\n", duration.count()); // TODO use gurobi Runtime model attr
        if (cb)
        {
            printf("Number of callbacks: %d\n", cb->numCallbacks);
            printf("Time in callbacks: %lf seconds\n", cb->callbackTime);
            printf("Number of lazy constraints generated: %d\n", cb->numLazyCuts);
            delete cb;
        }

        // will remain temporary for script run
        if(model.get(GRB_IntAttr_Status) == 3) // infeasible
          printf("qwerky567: model is infeasible\n");
        else {

//        double objval = model.get(GRB_DoubleAttr_ObjVal);
//        printf("qwerky567: Objective value: %lf (%e), time: %lf seconds, MIP gap: %.2lf%%, Bound: %lf\n", objval, objval, duration.count(), model.get(GRB_DoubleAttr_MIPGap)*100., model.get(GRB_DoubleAttr_ObjBound));

        // will remain temporary for script run
        //if(model.get(GRB_IntAttr_SolCount) > 0)
        //  printf("qwerky567: sol: L = %.4lf, U = %.4lf\n", x[g->nr_nodes][0].get(GRB_DoubleAttr_X), x[g->nr_nodes][1].get(GRB_DoubleAttr_X));
        //else
        //  printf("qwerly567: sol: no incumbent solution found!\n");
        }

        if (need_solution) {
            vector<int> sol;
            translate_solution(x, sol, g->nr_nodes);
            printf_solution(sol, "districting.out");
        }

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
