// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <cmath>

using namespace std;

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

// obsolete for assymetric formulation where only i < j
// only for j < i
/*inline int ind(const int i, const int j, const int n)
{
  return j + (2*n - i - 1)*i/2;
}*/

//for debug
void print_GRBLinExpr(GRBLinExpr& expr)
{
  for(unsigned int i = 0; i < expr.size(); ++i)
  {
    if(expr.getCoeff(i) != 0)
      printf("+ %.1lf %s", expr.getCoeff(i), expr.getVar(i).get(GRB_StringAttr_VarName).c_str());
  }
  printf("\n");
//  printf(" <= %.1lf", expr.getConstant());
}


class LazyCallback: public GRBCallback
{
  private:
    double* x;
    int* visited;
    int* aci;
    GRBVar* vars;
    int n;
    int nvars;
    graph* g;
  public:
    int numCallbacks;
    double callbackTime;
    LazyCallback(GRBVar* v, int n_, graph* g_) : n(n_), numCallbacks(0), callbackTime(0.)
    {
      nvars = n_*(n_+1)/2; x = new double[nvars];
      vars = v;
      visited = new int[n_];
      aci = new int[n_]; // who cares about memory, maximum time/space trade-off
      g = g_;
//      x = new double*[nrd_];
    }
#define mydel(x) { if(x) delete [] x; }
    virtual ~LazyCallback() {
      mydel(x);
      mydel(visited);
      mydel(aci);
    }
  protected:
    void callback () {
      try {
        if (where == GRB_CB_MIPSOL) {
//          fprintf(stderr, "Entered CALLBACK!\n");
          numCallbacks++;
          auto start = chrono::steady_clock::now();
          // MIP solution callback
//          for(int i = 0; i < nrd; ++i)
            x = getSolution(vars, nvars);
//          for(int i = 0; i < n; ++i)
//            sol[i] = (x[ind(i,i,n)] > 0.5);
/*          for(int i = 0; i < n; ++i)
            for(int j = 0; j <= i; ++j)
              printf("x%d_%d (%d) = %d\n", j, i, ind(j,i,n), (x[ind(j,i,n)]>0.5)?1:0);*/
/*          for(int i = 0; i < n; ++i)
            if(x[ind(i,i,n)] > 0.5)
            {
               printf("Cluster head %d:", i);
               for(int j = 0; j < i; ++j) // clusterhead is always maximum!
               {
               if(x[ind(j,i,n)] > 0.5)
               {
                  printf(" %u", j);
                }
             }
             printf("\n");
           }
*/

          // main work done here
          vector<int> s;
          // check constraint 7
          // run DFS from ii checking al j < i xji
          for(int i = 0; i < n; ++i)
          {
            if(x[ind(i,i,n)] > 0.5) // i is a cluster head
            {
              // run DFS on C_i and mark vertices
              for(int k = 0; k < n; ++k)
                visited[k] = 0;
              visited[i] = 1;
              s.clear(); s.push_back(i);
              while(!s.empty())
              {
                int cur = s.back(); s.pop_back();
                visited[cur] = 1;
                auto nb = g->nb(cur);
                for(auto it = nb.begin(); it != nb.end(); ++it)
                {
                  if(*it < i && x[ind(*it, i, n)] > 0.5) // in C_i
                  {
                    if(!visited[*it])
                      s.push_back(*it);
                  }
                }
              }
              for(int k = 0; k < i; ++k)
              {
                if(x[ind(k,i, n)] > 0.5 && !visited[k])
                {
//                  printf("Cluster head %d not connected to node %d\n", i, k);
                  // add k,i-separator now
                  // mark A(C_i), run DFS from i
                  for(int j = 0; j < n; ++j)
                  {
                    visited[j] = 0;
                    aci[j] = 0;
                  }
                  s.clear(); s.push_back(i);
                  while(!s.empty())
                  {
                    int cur = s.back(); s.pop_back();
                    visited[cur] = 1;
                    auto nb = g->nb(cur);
                    for(auto it = nb.begin(); it != nb.end(); ++it)
                    {
                      if(!visited[*it])
                      {
                        if(*it < i && x[ind(*it, i, n)] > 0.5) // *it in C_i
                          s.push_back(*it);
                        else
                          aci[*it] = 1;
                      }
                    }
                  }
/*                  fprintf(stderr,"Done marking A(C_i)\nA(C_i): ");
                  for(int j = 0; j < n; ++j)
                    if(aci[j])
                      printf("%d ", j);
                  printf("\n");*/
                  // run DFS from k to compute R_i \cap A(C_i)
                  s.clear();
                  GRBLinExpr expr = 0;
                  for(int j = 0; j < n; ++j)
                    visited[j] = 0;
                  s.push_back(k); visited[k] = 1;
                  while(!s.empty())
                  {
                    int cur = s.back(); s.pop_back();
                    auto nb = g->nb(cur);
                    for(auto it = nb.begin(); it != nb.end(); ++it)
                    {
                      if(!visited[*it])
                      {
                        visited[*it] = 1;
                        if(aci[*it])
                        {
                          if(*it < i)
                          {
                            //printf("+ x%d_%d (%d)", *it, i, ind(*it, i, n));
                            //printf("+ C%d", ind(*it, i, n));
                            expr += vars[ind(*it, i, n)];
                          }
                        } else {
                          s.push_back(*it);
                        }
                      }
                    }
                  }
                  //printf(" >= x%d_%d (%d)\n", k, i, ind(k,i,n));
                  //printf(" - C%d >= 0\n", ind(k,i,n));
                  //printf("Actual expression (>= 0): ");
                  expr -= vars[ind(k,i,n)];
                  //print_GRBLinExpr(expr);
                  addLazy(expr >= 0);
                  goto myend; //FIXME but thats's easier, no?
                }
              }
            }
          }

myend: ;
          // record callback time
          chrono::duration<double> d = chrono::steady_clock::now() - start;
          callbackTime += d.count();
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
    }
};

int main(int argc, char *argv[])
{
  ios::sync_with_stdio(1);
  if(argc < 3)
  {
    printf("Usage: %s <dimacs input> <data input>\n", argv[0]);
    return 0;
  }
  try {
    graph* g = from_dimacs(argv[1]);
    if(!g)
    {
      printf("Fail to load graph %s\n", argv[1]);
      return 1;
    }
    if(!g->is_connected())
    {
      printf("Problem is infeasible (not connected!)\n");
      return 0;
    }
    int n = g->nr_nodes;

    FILE* f = fopen(argv[2], "r");
    if(!f)
    {
      printf("Failed to open data file %s\n", argv[2]);
      return 1;
    }

    // read L, U
    int L = 0, U = 100000;
    // read k
    int k = 1; // number of district
    fscanf(f, "%d %d %d ", &k, &L, &U);

    // read p[i]
    // skip r b for compatibility with old gerry
    vector<int> p(n, 0); // or 1?
    for(int i = 0; i < n; ++i)
    {
      int tmp, tmp2, tmp3;
      fscanf(f, "%d %d %d ", &tmp, &tmp2, &tmp3); // ignore tmp2 tmp3, they are red and blue #
      p[i] = tmp;
    }
    // read dij
    vector<vector<int> > d(n, vector<int>(n, 1)); //TODO data with d_ij

    fclose(f);
    printf("starting gurobi. k = %d, L = %d, U = %d\n", k, L, U);

    GRBEnv env = GRBEnv();

    GRBModel model = GRBModel(env);
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    // Create variables
    int nr_vars = n*(n+1)/2; // store only i <= j
    GRBVar* x = model.addVars(nr_vars, GRB_BINARY);
    model.update();

    // Set objective: minimize sum x_ij
    GRBLinExpr expr = 0;
    for(int i = 0; i < n; ++i)
      for(int j = i; j < n; ++j)
        expr += d[i][j] * p[i] * x[ind(i,j,n)];

    model.setObjective(expr, GRB_MINIMIZE);

    // add constraints (2)
    for(int i = 0; i < n; ++i)
    {
      GRBLinExpr constr = 0;
      for(int j = i; j < n; ++j)
        constr += x[ind(i,j,n)];
      model.addConstr(constr == 1);
    }

    // add contraints (5)
    for(int i = 0; i < n; ++i)
    {
      for(int j = 0; j < i; ++j)
        model.addConstr(x[ind(j,i,n)] <= x[ind(i,i,n)]);
    }

    // add constraint 3
    expr = 0;
    for(int j = 0; j < n; ++j)
      expr += x[ind(j,j,n)];
    model.addConstr(expr == k);

    // add constraint 4
    for(int j = 0; j < n; ++j)
    {
      GRBLinExpr constr = 0;
      for(int i = 0; i < j; ++i)
        constr += p[i]*x[ind(i,j,n)];
      // add for j
      model.addConstr(constr + (p[j] - U)*x[ind(j,j,n)] <= 0); // U
      model.addConstr(constr + (p[j] - L)*x[ind(j,j,n)] >= 0); // L
    }

    // Optimize model
    LazyCallback cb = LazyCallback(x, n, g);
    model.setCallback(&cb);
    model.write("out.lp");
    model.optimize();

    // translate solution
    for(int i = 0; i < n; ++i)
      if(x[ind(i,i,n)].get(GRB_DoubleAttr_X) > 0.5)
      {
        printf("Cluster head %d:", i);
        for(int j = 0; j < i; ++j) // clusterhead is always maximum!
        {
          if(x[ind(j,i,n)].get(GRB_DoubleAttr_X) > 0.5)
          {
            printf(" %u", j);
          }
        }
        printf("\n");
      }
//      for(int j = i; j < n; ++j)
//        printf("x[%d][%d] = %lf\n", i, j, x[ind(i,j,n)].get(GRB_DoubleAttr_X));

    cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    cout << "Number of callbacks = " << cb.numCallbacks << endl;
    cout << "Total time in callbacks = " << cb.callbackTime << " secs " << endl;

    delete g;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
