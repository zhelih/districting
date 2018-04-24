#include "gurobi_c++.h"
#include "KGraph.h"
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

string itos(int i);
string statusNumtoString(int num);
vector<vector<long>> solveMST(KGraph &g1, KGraph &g2);
vector<vector<long>> solveFM2(KGraph &g1, KGraph &g2, KGraph &d, long &K, long &L, long &U);
vector<vector<long>> solveFM3(KGraph &g1, KGraph &g2, KGraph &d, long &K, long &L, long &U);
