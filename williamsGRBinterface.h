#include "gurobi_c++.h"
#include "KGraph.h"
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

string itos(int i);
string statusNumtoString(int num);
vector<vector<long>> solveMST(KGraph &g1, KGraph &g2);
