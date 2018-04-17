// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include <cstdio>
#include <chrono>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_set> // for variable reduction
#include <queue> // for Dijkstra
#include <climits> // for INT_MAX

using namespace std;

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef max
#define max(a,b) (((a)<(b))?(b):(a))
#endif

// obsolete for assymetric formulation where only i < j
// only for j < i
inline int ind(const int i, const int j, const int n)
{
  return j + (2*n - i - 1)*i/2;
}

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

// magic going on here
// required by unordered_set
namespace std {
template <> struct hash<std::pair<int, int>> {
    inline size_t operator()(const std::pair<int, int> &v) const {
        std::hash<int> int_hasher;
        return int_hasher(v.first) ^ int_hasher(v.second);
    }
};

}

// @return number of variables removed
int populate_zeros(unordered_set<pair<int,int>> &x_zeros, graph* g, const vector<int> &p, int U);

class LazyCallback: public GRBCallback
{
  private:
    int** x;
    int* visited;
    int* aci;
    GRBVar*** vars;
    int n;
    graph* g;
  public:
    int numCallbacks;
    double callbackTime;
    LazyCallback(GRBVar*** v, int n_, graph* g_) : n(n_), numCallbacks(0), callbackTime(0.)
    {
      x = new int*[n];
      for(int i = 0; i < n; ++i)
        x[i] = new int[n];
      vars = v; // pointer copy
      visited = new int[n_];
      aci = new int[n_]; // who cares about memory, maximum time/space trade-off
      g = g_; // shallow copy!
    }
#define mydel(x) { if(x) delete [] x; }
    virtual ~LazyCallback() {
      for(int i = 0; i < n; ++i)
        mydel(x[i]);
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
          for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j)
              x[i][j] = (vars[i][j] && getSolution(vars[i][j][0]) > 0.5)?1:0;

          // main work done here
          vector<int> s;
          // check constraint 7
          // run DFS from ii
          for(int i = 0; i < n; ++i)
          {
            if(x[i][i]) // i is a cluster head
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
                auto nbs = g->nb(cur);
                for(int nb: nbs)
                {
                  if(x[nb][i]) // in C_i
                  {
                    if(!visited[nb])
                      s.push_back(nb);
                  }
                }
              }
              for(int k = 0; k < n; ++k)
              {
                if(x[k][i] && !visited[k])
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
                    auto nbs = g->nb(cur);
                    for(int nb: nbs)
                    {
                      if(!visited[nb])
                      {
                        if(x[nb][i]) // *it in C_i
                          s.push_back(nb);
                        else
                          aci[nb] = 1;
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
                    auto nbs = g->nb(cur);
                    for(int nb: nbs)
                    {
                      if(!visited[nb])
                      {
                        visited[nb] = 1;
                        if(aci[nb])
                        {
                          {
                            //printf("+ x%d_%d (%d)", *it, i, ind(*it, i, n));
                            //printf("+ C%d", ind(*it, i, n));
                            expr += vars[nb][i][0];
                          }
                        } else {
                          s.push_back(nb);
                        }
                      }
                    }
                  }
                  //printf(" >= x%d_%d (%d)\n", k, i, ind(k,i,n));
                  //printf(" - C%d >= 0\n", ind(k,i,n));
                  //printf("Actual expression (>= 0): ");
                  expr -= vars[k][i][0];
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
    printf("Usage: %s [ <dimacs input> <data input> | grid <n> ]\n", argv[0]);
    return 0;
  }
  try {
    graph* g = 0;
    int L = 0, U = 100000;
    int k = 1; // number of district
    vector<int> p;
    if(!strcmp(argv[1], "grid"))
    {
      int n = atoi(argv[2]);
      g = from_grid(n);
      L = n; U = n; k = n;
      p.resize(n*n);
      for(int i = 0; i < n*n; ++i)
        p[i] = 1;

    } else {
      g = from_dimacs(argv[1]);
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

      // skip first two lines
      char buf[1023];
      fgets(buf, sizeof(buf), f);
      fgets(buf, sizeof(buf), f);

      // read L, U and read k
      fscanf(f, "%d %d %d ", &k, &L, &U);

      // read p[i]
      p.resize(n);
      for(int i = 0; i < n; ++i)
      {
        int tmp, tmp2;
        fscanf(f, "%d %d ", &tmp, &tmp2);
        p[tmp] = tmp2;
      }
      fclose(f);

    }
    int n = g->nr_nodes;
    vector<vector<float> > d(n, vector<float>(n, 1)); //TODO data with d_ij
    // Euclidian distances for a grid
/*      int n_ = sqrt(n);
      for(int i1 = 0; i1 < n_; ++i1)
        for(int j1 = 0; j1 < n_; ++j1)
          for(int i2 = 0; i2 < n_; ++i2)
            for(int j2 = 0; j2 < n_; ++j2) {
              int delta_i = abs(i1-i2);
              int delta_j = abs(j1-j2);
              float d_ = sqrt((float)(delta_i*delta_i + delta_j*delta_j)); // euclidian
              float d2_ = delta_i + delta_j; // hop (manhattan)
              float d3_ = max(delta_i, delta_j); // chebyshev (l_infinity)
              d[i1+n_*j1][i2+n_*j2] = d2_;
            }
*/
    g->floyd_warshall(d);

    // reduce number of variables (safe)
    unordered_set<pair<int,int>> x_zeros;
    int rem = populate_zeros(x_zeros, g, p, U);
    printf("Preprocessing : removed %d variables (U)\n", rem);
    printf("starting gurobi. k = %d, L = %d, U = %d\n", k, L, U);

    GRBEnv env = GRBEnv();

    GRBModel model = GRBModel(env);
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
    // Create variables ( wut )
    GRBVar*** x = new GRBVar**[n];
    for(int i = 0; i < n; ++i)
      x[i] = new GRBVar*[n];

    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(x_zeros.count(make_pair(i,j)) > 0)
          x[i][j] = NULL;
        else
          x[i][j] = model.addVars(1, GRB_BINARY);
    model.update();

    // Set objective: minimize sum x_ij
    GRBLinExpr expr = 0;
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(x[i][j])
          expr += d[i][j] * p[i] * x[i][j][0];

    model.setObjective(expr, GRB_MINIMIZE);

    // add constraints (2)
    for(int i = 0; i < n; ++i)
    {
      GRBLinExpr constr = 0;
      for(int j = 0; j < n; ++j)
        if(x[i][j])
          constr += x[i][j][0];
      model.addConstr(constr == 1);
    }

    // add contraints (5)
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(i != j && x[i][j]) model.addConstr(x[j][i][0] <= x[i][i][0]);

    // add constraint 3
    expr = 0;
    for(int j = 0; j < n; ++j)
      expr += x[j][j][0];
    model.addConstr(expr == k);

    // add constraint 4
    for(int j = 0; j < n; ++j)
    {
      GRBLinExpr constr = 0;
      for(int i = 0; i < n; ++i)
        if(i != j && x[i][j]) constr += p[i]*x[i][j][0];
      // add for j
      model.addConstr(constr + (p[j] - U)*x[j][j][0] <= 0); // U
      model.addConstr(constr + (p[j] - L)*x[j][j][0] >= 0); // L
    }

    // Optimize model
    LazyCallback cb = LazyCallback(x, n, g);
    model.setCallback(&cb);
    model.write("out.lp");
    model.optimize();

    // translate solution
    for(int i = 0; i < n; ++i)
      if(x[i][i][0].get(GRB_DoubleAttr_X) > 0.5)
      {
        printf("Cluster head %d:", i);
        for(int j = 0; j < n; ++j)
        {
          if(j != i && x[j][i] && x[j][i][0].get(GRB_DoubleAttr_X) > 0.5)
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

    // write to file for an output
    vector<vector<int> > clusters(n);;
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(x[i][j] && x[i][j][0].get(GRB_DoubleAttr_X) > 0.5)
          clusters[j].push_back(i);

/*    printf("Districts:\n");
    for(i = 0; i < nr_districts; ++i)
    {
      printf("District %d: (y[i] = %d)", i, ((y[i].get(GRB_DoubleAttr_X)>0.5)?1:0));
      for(auto it = clusters[i].begin(); it != clusters[i].end(); ++it)
        printf(" %d", *it);
      printf("\n");
    }
*/
    FILE* f = fopen("districting.out", "w");
    int nr = 0;
    for(int i = 0; i < n; ++i)
    {
      for(int node : clusters[i])
        fprintf(f, "%d %d\n", node, nr);
      nr += (clusters[i].size() > 0);
    }

    fclose(f);


    delete g;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}

// @return number of variables removed
int populate_zeros(unordered_set<pair<int,int>> &x_zeros, graph* g, const vector<int> &p, int U) {
  // run Dijsktra with d(s) = p[s] and w[u,v] = p[v]
  // scan all distances, if legnth > U then add to x_zeros
  // TODO optimize time? Is this really takes long to spend time on optimizing?
  for(int s = 0; s < g->nr_nodes; ++s) { // select s as source
    // typical Dijsktra with heap here
    priority_queue< pair<int,int>, vector <pair<int,int>> , greater<pair<int,int>> > pq;
    vector<int> dist(g->nr_nodes, INT_MAX); //FIXME long?
    pq.push(make_pair(0, s)); // copy constructor?
    dist[s] = p[s]; // NB: not zero here!
    while (!pq.empty())
    {
        int u = pq.top().second;
        pq.pop();
        for(int nb: g->nb(u)) {
          int weight = p[nb];
          if(dist[nb] > dist[u] + weight) {
            dist[nb] = dist[u] + weight;
            pq.push(make_pair(dist[nb], nb));
          }
        }
    }

    for(int i = 0; i < g->nr_nodes; ++i) {
      if(dist[i] > U)
        x_zeros.insert(make_pair(s, i));
    }
  }
/*  for(pair<int,int> p : x_zeros) {
    printf("(%d,%d) ", p.first, p.second);
  }
  printf("\n");*/
  return x_zeros.size();
}
