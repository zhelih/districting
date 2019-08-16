#ifndef _IO_H
#define _IO_H
#include "graph.h"
#include "models.h"
#include <vector>
#include <string>

using namespace std;

struct run_params
{
  string dimacs_file;
  string distance_file;
  string population_file;
  char state[3];
  int L;
  int U;
  int k;
  string model;
  string ralg_hot_start;
  string output;
};

run_params read_config(const char* fname, const char* state, const char* ralg_hot_start);

int read_input_data(const char* dimacs_fname, const char* distance_fname, const char* population_fname, // INPUTS
                     graph* &g, vector<vector<int> >& dist, vector<int>& population); // OUTPUTS
// construct districts from hess variables
void translate_solution(hess_params& p, vector<int>& sol, int n);
// prints the solution <node> <district>
void printf_solution(const vector<int>& sol, const char* fname=NULL);
void calculate_UL(const vector<int> population, int k, int* L, int* U);
int read_auto_int(const char*, int);
//read ralg initial point from file [fname] to [x0]
void read_ralg_hot_start(const char* fname, double* x0, int dim);
// dump result to "ralg_hot_start.txt"
void dump_ralg_hot_start(double* res, int dim);

#endif
