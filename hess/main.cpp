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
#include <string>
#include "ralg/ralg.h"

#ifndef sign
#define sign(x) (((x)>0)?1:((x)==0)?0:(-1))
#endif

using namespace std;
//extern const char* gitversion;

#define DO_BATCH_OUTPUT

int main(int argc, char *argv[])
{
    ios::sync_with_stdio(1);
    //printf("Districting, build %s\n", gitversion);
    if (argc < 8) {
        printf("Usage: %s <dimacs> <distance> <population> <L|auto> <U|auto> <k> <model> [ralg hot start]\n\
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
    bool ralg_hot_start = argc > 8;
    char* ralg_hot_start_fname = ralg_hot_start?argv[8]:nullptr;

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
    printf("Model input: total population = %ld\n", total_pop);

    string arg_model = argv[7];

    //apply the merging preprocess and get the clusters
    vector<vector<int>> clusters;
    if (arg_model == "scf" || arg_model == "mcf0" || arg_model == "mcf1" || arg_model == "mcf2" || arg_model == "cut1" || arg_model == "cut2")
    {
        printf("Preprocessing the graph...\n");
        vector<int> new_population;
        new_population = population;
        clusters = preprocess(g, new_population, L, U, population);
    }

    //apply Lagrangian algorithm
    vector<vector<double>> w(g->nr_nodes, vector<double>(g->nr_nodes)); // this is the weight matrix in the objective function

    for (int i = 0; i < g->nr_nodes; i++)
        for (int j = 0; j < g->nr_nodes; j++)
            w[i][j] = get_objective_coefficient(dist, population, i, j);

    vector<vector<bool>> F_0(g->nr_nodes, vector<bool>(g->nr_nodes, false)); // define matrix F_0
    vector<vector<bool>> F_1(g->nr_nodes, vector<bool>(g->nr_nodes, false)); // define matrix F_1

    vector<bool> S(g->nr_nodes);

    vector<vector<double>> w_hat(g->nr_nodes, vector<double>(g->nr_nodes));
    vector<double> W(g->nr_nodes, 0);


    // fix clusters to singletones for now
    clusters.clear();
    clusters.resize(g->nr_nodes);
    for (int i = 0; i < g->nr_nodes; ++i)
        clusters[i].push_back(i);

    auto cb_grad_func = [g, L, U, k, &population, &w, &S, &F_0, &F_1, &W, &w_hat, &clusters](const double* multipliers, double& f_val, double* grad) {
        f_val = 0;

        // calculate here the gradient and obj value
        solveInnerProblem(g, multipliers, F_0, F_1, L, U, k, clusters, population, w, w_hat, W, S, grad, f_val);

        // revert the grads if needed to maintain positive l,u
        for (int i = g->nr_nodes; i < 3 * g->nr_nodes; ++i)
            grad[i] = sign(multipliers[i])*grad[i];
        return true;
    };

    double UB = INFINITY; // = 1439.05;
    vector<long> frequency(g->nr_nodes, 0);
    double bestLB = -INFINITY;
    vector<bool> bestCenters(g->nr_nodes, false);
    vector<double> bestW(g->nr_nodes, 0);
    vector<vector<double>> bestw_hat(g->nr_nodes, vector<double>(g->nr_nodes));
    vector<double> bestMultipliers(3*g->nr_nodes, 0);

    auto cb_grad_func2 = [g, L, U, k, &population, &w, &W, &w_hat, &S, &F_0, &F_1, &clusters, &UB, &frequency, &bestLB, &bestCenters, &bestW, &bestw_hat](const double* multipliers, double& f_val, double* grad) {
        f_val = 0;

        // calculate here the gradient and obj value
        eugene_inner(g, multipliers, L, U, k, population, w, w_hat, W, grad, f_val, S, F_0, F_1);
        //lagrangianBasedSafeFixing(F_0, F_1, clusters, W, S, f_val, UB);
        if (f_val > bestLB)
        {
            bestLB = f_val;
            bestCenters = S;
            bestW = W;
        }

        for (int i = 0; i < g->nr_nodes; ++i)
        {
            if (S[i]) frequency[i]++;
        }

        for (int i = 0; i < g->nr_nodes; ++i) bestw_hat[i] = w_hat[i];

        return true;
    };

    auto start = chrono::steady_clock::now();

    // run ralg
    int dim = 3 * g->nr_nodes;
    double * x0 = new double[dim];
    // try to load hot start if any
    if(ralg_hot_start)
      read_ralg_hot_start(ralg_hot_start_fname, x0, dim);
    else
      for (int i = 0; i < dim; ++i)
        x0[i] = 1.; // whatever
    double* res = new double[dim];
    ralg_options opt = defaultOptions; opt.output_iter = 1;
    double LB = ralg(&opt, cb_grad_func2, dim, x0, res, RALG_MAX); // lower bound from lagrangian

    // dump result to "ralg_hot_start.txt"
    dump_ralg_hot_start(res, dim);

    delete[] x0;
    delete[] res;

    /*sort(W.begin(), W.end());
    cerr << "After ralg, the final W values are: ";
    for (int i = 0; i < g->nr_nodes; ++i) cerr << W[i] << " ";
    cerr << endl;*/

    /*cerr << "And, the set of centers is (k=" << k << ") : ";
    for (int i = 0; i < g->nr_nodes; ++i)
        if (S[i])
            cerr << i << " ";
    cerr << endl;*/

    //cerr << "total number of variables = " << g->nr_nodes*g->nr_nodes << " =? "<< numFixedZero + numFixedOne + numUnfixed << endl;

    //FIXME floods the output
    /*for (int i = 0; i < g->nr_nodes; i++)
        cout << "vertex " << i << " , " << S[i] << endl;*/

    try {
        // initialize environment and create an empty model
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        bool need_solution = true;

        // get incumbent solution using centers from lagrangian
        GRBVar** x = 0;
        if (arg_model != "ul1" && arg_model != "ul2")
            x = build_hess(&model, g, dist, population, L, U, k);

        // push GUROBI to branch over clusterheads
        for (int i = 0; i < g->nr_nodes; ++i)
            x[i][i].set(GRB_IntAttr_BranchPriority, 1);

        HessCallback* cb = 0;

        //if (arg_model != "hess" && arg_model != "ul1" && arg_model != "ul2")
        //	strengthen_hess(&model, x, g, clusters);

        if (arg_model == "scf")
            build_scf(&model, x, g);
        else if (arg_model == "mcf0")
            build_mcf0(&model, x, g);
        else if (arg_model == "mcf1")
            build_mcf1(&model, x, g);
        else if (arg_model == "mcf2")
            build_mcf2(&model, x, g);
        else if (arg_model == "cut1")
            cb = build_cut1(&model, x, g);
        else if (arg_model == "cut2")
            cb = build_cut2(&model, x, g);
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


        //  solve restricted IP to get an upper bound 
        //		first, fix the non-selected center vars x[j][j] to zero.
        //for (int j = 0; j < g->nr_nodes; ++j) if (!S[j]) x[j][j].set(GRB_DoubleAttr_UB, 0);

        // select k smallest
        vector<int> frequency_indices(frequency.size());
        for (size_t i = 0; i < W.size(); ++i)
            frequency_indices[i] = i;
        sort(frequency_indices.begin(), frequency_indices.end(), [&frequency](int i1, int i2) { return frequency[i1] < frequency[i2]; });

        cerr << "frequency values: ";
        for (int i = 0; i < g->nr_nodes; ++i)
        {
            int v = frequency_indices[i];
            cerr << frequency[v] << " ";
        }
        cerr << endl;

        /*for (int i = 0; i < g->nr_nodes-2*k; ++i)
        {
            int v = frequency_indices[i];
            x[v][v].set(GRB_DoubleAttr_UB, 0);
        }*/

        cerr << "frequency values for best centers: ";
        for (int i = 0; i < g->nr_nodes; ++i)
        {
            if (!bestCenters[i]) x[i][i].set(GRB_DoubleAttr_UB, 0);
            else cerr << frequency[i] << " ";
        }
        cerr << endl;

        /*for (int i = 0; i < g->nr_nodes; ++i)
        {
            if(frequency[i]==0) x[i][i].set(GRB_DoubleAttr_UB, 0);
        }*/

        model.optimize();
        double UB = model.get(GRB_DoubleAttr_ObjVal);
        cerr << "UB from restricted IP = " << UB << endl;

        // unfix the vars that had been fixed (in the restricted problem)
        for (int j = 0; j < g->nr_nodes; ++j)
        {
            x[j][j].set(GRB_DoubleAttr_UB, 1);
            x[j][j].set(GRB_DoubleAttr_LB, 0);
        }

        // figure out which x[i][j] vars can be fixed to zero/one and fix them.
        lagrangianBasedSafeFixing(F_0, F_1, clusters, bestW, bestCenters, bestLB, UB, bestw_hat);
        for (int i = 0; i < g->nr_nodes; ++i)
        {
            for (int j = 0; j < g->nr_nodes; ++j)
            {
                if (F_0[i][j]) x[i][j].set(GRB_DoubleAttr_UB, 0);
                if (F_1[i][j]) x[i][j].set(GRB_DoubleAttr_LB, 1);
            }
        }

        //TODO change user-interactive?
        model.set(GRB_DoubleParam_TimeLimit, 3600.); // 1 hour
        model.set(GRB_IntParam_Threads, 10);
        model.set(GRB_DoubleParam_NodefileStart, 10); // 10 GB
        model.set(GRB_IntParam_Method, 3);
        model.set(GRB_DoubleParam_MIPGap, 0);

        //optimize the model

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

#ifdef DO_BATCH_OUTPUT

        printf("qwerky567: %s, %d, %d, %d, %d, %.2lf", dimacs_fname, k, g->nr_nodes, L, U, duration.count());

        // output overtly
        int max_pv = population[0];
        for (int pv : population)
            max_pv = max(max_pv, pv);


        printf(",%.2lf", static_cast<double>(max_pv) / static_cast<double>(U));

        // will remain temporary for script run
        if (model.get(GRB_IntAttr_Status) == 3) // infeasible
            printf(",infeasible,,");
        else {

            double objval = model.get(GRB_DoubleAttr_ObjVal);
            double mipgap = model.get(GRB_DoubleAttr_MIPGap)*100.;
            double objbound = model.get(GRB_DoubleAttr_ObjBound);


            // no incumbent solution was found, these values do no make sense
            if (model.get(GRB_IntAttr_SolCount) == 0)
            {
                mipgap = 100.;
                objval = 0.;
            }

            printf(", %.2lf, %.2lf, %.2lf", objval, mipgap, objbound);
        }

        long nodecount = static_cast<long>(model.get(GRB_DoubleAttr_NodeCount));

        int num_callbacks = 0;
        double time_callbacks = 0.;
        int num_lazy = 0;
        if (cb) {
            num_callbacks = cb->numCallbacks;
            time_callbacks = cb->callbackTime;
            num_lazy = cb->numLazyCuts;
        }
        // state, k n l u time obj mipgap objbound nodes callback x3
        printf(", %ld, %d, %.2lf, %d\n", nodecount, num_callbacks, time_callbacks, num_lazy);

        // will remain temporary for script run
        //if(model.get(GRB_IntAttr_SolCount) > 0)
        // printf("qwerky567: sol: L = %.4lf, U = %.4lf\n", x[g->nr_nodes][0].get(GRB_DoubleAttr_X), x[g->nr_nodes][1].get(GRB_DoubleAttr_X));
        //else
        // printf("qwerly567: sol: no incumbent solution found!\n");
#endif

        if (need_solution && model.get(GRB_IntAttr_Status) != 3) {
            vector<int> sol;
            translate_solution(x, sol, g->nr_nodes);
            string fn = string(dimacs_fname);
            string soln_fn = fn.substr(0, 2) + "_" + arg_model + ".sol";
            int len = soln_fn.length() + 1;
            char* char_array = new char[len];
            strcpy(char_array, soln_fn.c_str());
            printf_solution(sol, char_array);
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
