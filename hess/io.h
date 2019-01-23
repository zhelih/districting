#ifndef _IO_H
#define _IO_H
#include "graph.h"
#include "gurobi_c++.h"
#include <vector>

using namespace std;

int read_input_data(const char* dimacs_fname, const char* distance_fname, const char* population_fname, // INPUTS
                     graph* &g, vector<vector<int> >& dist, vector<int>& population); // OUTPUTS
// construct districts from hess variables
void translate_solution(GRBVar** x, vector<int>& sol, int n);
// prints the solution <node> <district>
void printf_solution(const vector<int>& sol, const char* fname=NULL);
void calculate_UL(const vector<int> population, int k, int* L, int* U);
int read_auto_int(const char*, int);
#endif
