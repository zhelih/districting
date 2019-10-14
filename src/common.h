#ifndef __COMMON_H
#define __COMMON_H

#include <vector>
#include <unordered_map>
#include "gurobi_c++.h"

struct hess_params
{
  GRBVar* x;
  std::vector<std::vector<bool>> F0;
  std::vector<std::vector<bool>> F1;
  std::unordered_map<int, int> h;
  int n;
  int infty;
};

//hack
#define IS_X(i,j) (!p.F0[i][j] && !p.F1[i][j])
#define X_V(i,j) (p.x[p.h[p.n*i+j]])
#define X(i,j) (p.F0[i][j]?GRBLinExpr(0.):(p.F1[i][j]?GRBLinExpr(1.):GRBLinExpr(X_V(i,j))))

struct run_params
{
  std::string dimacs_file;
  std::string distance_file;
  std::string population_file;
  char state[3];
  int L;
  int U;
  int k;
  std::string model;
  std::string ralg_hot_start;
  FILE* output;
};

struct pair_hash {
  inline std::size_t operator()(const std::pair<int, int> & v) const {
    return v.first * 31 + v.second;
  }
};

#define mymin(a,b) (((a)<(b))?(a):(b))
#define mymax(a,b) (((a)>(b))?(a):(b))
#define myabs(x) (((x)>0)?(x):(-(x)))

#define MYINFINITY 1e20

#endif
