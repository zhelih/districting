// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include "io.h"
#include "models.h"
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
using namespace std;

extern const char* gitversion;


int main(int argc, char *argv[])
{
  ios::sync_with_stdio(1);
  printf("Districting, build %s\n", gitversion);
  if(argc < 3) {
    printf("Usage: %s <dimacs> <distance> <population> <L> <U> <k>\n", argv[0]);
    return 0;
  }
  // parse command line arguments
  char* dimacs_fname = argv[1];
  char* distance_fname = argv[2];
  char* population_fname = argv[3];
  int L = std::stoi(argv[4]);
  int U = std::stoi(argv[5]);
  int k = std::stoi(argv[6]);
  printf("Model input: L = %d, U = %d, k = %d\n", L, U, k);

  // read inputs
  graph* g = 0;
  vector<vector<int> > dist;
  vector<int> population;
  if(read_input_data(dimacs_fname, distance_fname, population_fname, g, dist, population))
    return 1; // failure

  // check connectivity
  if(!g->is_connected())
  {
    printf("Problem is infeasible (not connected!)\n");
    return 1;
  }

  if(g->nr_nodes <= 0)
  {
    fprintf(stderr, "empty graph\n");
    return 1;
  }

  if(dist.size() != g->nr_nodes || population.size() != g->nr_nodes)
  {
    fprintf(stderr, "dist/population size != n, expected %d\n", g->nr_nodes);
    return 1;
  }

  try {
    // initialize environment and create an empty model
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // build hess model
    GRBVar** x = build_hess(&model, g, dist, population, L, U, k);

    // add SCF contraints
    build_scf(&model, x, g);

    //optimize the model
    model.optimize();

    vector<int> sol;
    translate_solution(x, sol, g->nr_nodes);

    printf_solution(sol, "districting.out");
  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(const char* msg) {
    cout << "Exception with message : " << msg << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }


  delete g;
  return 0;
}
