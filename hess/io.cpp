#include "graph.h"
#include "io.h"
#include <vector>
#include <cstdio>
#include <cstring>
#include <cmath>
#include "gurobi_c++.h"

using namespace std;

int read_input_data(const char* dimacs_fname, const char* distance_fname, const char* population_fname, // INPUTS
                     graph* &g, vector<vector<int> >& dist, vector<int>& population) // OUTPUTS
{
    // read dimacs graph
    g = from_dimacs(dimacs_fname);
    if(!g) {
      fprintf(stderr, "Failed to read dimacs graph from %s\n", dimacs_fname);
      return 1;
    }

    // read distances (must be sorted)
    FILE* f = fopen(distance_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", distance_fname);
      return 1;
    }
    // file contains the first row as a header row, skip it
    // also skip the first element in each row (node id)
    char buf[50000]; //dummy
    fgets(buf, sizeof(buf), f); // skip first line
    if(strlen(buf) >= sizeof(buf)-5)
      printf("WARNING: Possible buffer overlow!\n");
    dist.resize(g->nr_nodes);
    for(int i = 0; i < g->nr_nodes; ++i)
    {
      dist[i].resize(g->nr_nodes);
      int d; fscanf(f, "%d,", &d); // skip first element
      for(int j = 0; j < g->nr_nodes; ++j) {
        fscanf(f, "%d,", &d);
        dist[i][j] = d;
      }
    }
    fclose(f);

    // read population file
    f = fopen(population_fname, "r");
    if(!f) {
      fprintf(stderr, "Failed to open %s\n", population_fname);
      return 1;
    }
    // skip first line about total population
    fgets(buf, sizeof(buf), f);
    population.resize(g->nr_nodes);
    for(int i = 0; i < g->nr_nodes; ++i) {
      int node, pop;
      fscanf(f, "%d %d ", &node, &pop);
      population[node] = pop;
    }
    fclose(f);
    return 0;
}

// construct districts from hess variables
void translate_solution(hess_params& p, vector<int>& sol, int n)
{
    // translate the solution
    sol.resize(n);

    vector<int> heads(n, 0);
    int cur = 1;
    // firstly assign district number for clusterheads
    for(int i = 0; i < n; ++i)
    {
      if(p.F0[i][i])
        continue;
      if(p.F1[i][i] || X_V(i,i).get(GRB_DoubleAttr_X) > 0.5)
        heads[i] = cur++;
    }
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
        if(p.F0[i][j]) continue; else if (p.F1[i][j] || X_V(i,j).get(GRB_DoubleAttr_X) > 0.5)
          sol[i] = heads[j];
}

// prints the solution <node> <district>
void printf_solution(const vector<int>& sol, const char* fname)
{
  // write solution to file
  if(sol.size() > 0)
  {
    FILE* f;
    if(fname)
      f = fopen(fname, "w");
    else
      f = stderr;

    for(int i = 0; i < sol.size(); ++i)
      fprintf(f, "%d %d\n", i, sol[i]);

    if(fname)
      fclose(f);
  }
}

void calculate_UL(const vector<int> population, int k, int* L, int* U)
{
  int total_pop = 0;
  for(int p : population)
    total_pop += p;
  if(*U == 0)
     *U = static_cast<int>(floor(static_cast<float>(total_pop) / static_cast<float>(k) * 1.005));
  if(*L == 0)
    *L = static_cast<int>(ceil(static_cast<float>(total_pop) / static_cast<float>(k) * 0.995));
}

int read_auto_int(const char* arg, int def)
{
  if(strcmp(arg, "auto"))
    return std::atoi(arg);
  return def;
}

//read ralg initial point from file [fname] to [x0]
void read_ralg_hot_start(const char* fname, double* x0, int dim)
{
  FILE* f = fopen(fname, "r");
  if(!f)
  {
    fprintf(stderr, "Failed to open %s.\n", fname);
    return;
  }
  for(int i = 0; i < dim; ++i)
  {
    double val = 0.;
    if(EOF == fscanf(f, "%lf ", &val))
    {
      fprintf(stderr, "Failure to complete reading ralg x0 from %s.\n", fname);
      break;
    }
    x0[i] = val;
  }
  fclose(f);
}

// dump result to "ralg_hot_start.txt"
void dump_ralg_hot_start(double* res, int dim)
{
  const char* outname = "ralg_hot_start.txt";
  FILE* f = fopen(outname, "w");
  if(!f)
  {
    fprintf(stderr, "Cannot open %s for dumping ralg result.\n", outname);
    return;
  }
  for(int i = 0; i < dim; ++i)
    fprintf(f, "%lf\n", res[i]);
  fclose(f);
}

const char* parse_param(const char* src, const char* prefix)
{
  // check if src starts with prefix
  int len_prefix = strlen(prefix);

  if(len_prefix > strlen(src) + 2)
    return nullptr;

  for(int i = 0; i < len_prefix; ++i)
    if(src[i] != prefix[i])
      return nullptr;

  if(src[len_prefix] != ' ')
    return nullptr;

  int k = len_prefix;
  while(src[++k] == ' ');
  return src+k;
}

run_params read_config(const char* fname, const char* state)
{
  FILE* f = fopen(fname, "r");
  if(!f)
  {
    fprintf(stderr, "Cannot open config %s.\n", fname);
    exit(1);
  }

  run_params rp;

  if(strlen(state) > 1)
    strncpy(rp.state, state, 2);

  char buf[1020];
  string database;
  string level;
  bool seen_state = false;

  auto check_database = [&database]() { if(!database.empty()) { fprintf(stderr, "Config error: got dimacs/distance/population with database.\n"); exit(1); } };

  while(fgets(buf, sizeof(buf), f) != NULL)
  {
    if(strlen(buf) < 2 || buf[0] == '#')
      continue;

    const char* v;
    if((v = parse_param(buf, "database")) != nullptr)
      database = v; // copy
    else if((v = parse_param(buf, "level")) != nullptr)
      level = v;
    else if((v = parse_param(buf, "state")) != nullptr)
    {
      if(state != nullptr && state != NULL && strlen(state) > 0)
      {
        fprintf(stderr, "Found state in config file, however passed %s, which do you want me to use?\n", state);
        exit(1);
      }
      strncpy(rp.state, v, 2);
      seen_state = true;
    }
    else if((v = parse_param(buf, "dimacs")) != nullptr)
    {
      check_database();
      rp.dimacs_file = v;
    }
    else if((v = parse_param(buf, "distance")) != nullptr)
    {
      check_database();
      rp.distance_file = v;
    }
    else if((v = parse_param(buf, "population")) != nullptr)
    {
      check_database();
      rp.population_file = v;
    }
    else if((v = parse_param(buf, "model")) != nullptr)
      rp.model = v;
    else if((v = parse_param(buf, "ralg_hot_start")) != nullptr)
      rp.ralg_hot_start = v;
    else if((v = parse_param(buf, "output")) != nullptr)
      rp.output = v;
    else if((v = parse_param(buf, "L")) != nullptr)
    {
      if(strncmp(v, "auto", 4) == 0)
        rp.L = 0;
      else
        rp.L = atoi(v);
    }
    else if((v = parse_param(buf, "U")) != nullptr)
    {
      if(strncmp(v, "auto", 4) == 0)
        rp.U = 0;
      else
        rp.U = atoi(v);
    }
    else if((v = parse_param(buf, "k")) != nullptr)
    {
      if(strncmp(v, "auto", 4) == 0)
        rp.k = 0;
      else
        rp.k = atoi(v);
    }
  }
  fclose(f);

  // clean newlines from fgets
  auto clean_nl = [](string& s) { while(!s.empty() && s.back() == '\n') s.pop_back(); };

  clean_nl(database);
  clean_nl(level);
  clean_nl(rp.dimacs_file);
  clean_nl(rp.population_file);
  clean_nl(rp.distance_file);
  clean_nl(rp.model);
  clean_nl(rp.ralg_hot_start);
  clean_nl(rp.output);
  rp.state[2] = '\0';

  if(database.empty() && (rp.dimacs_file.empty() || rp.population_file.empty() || rp.distance_file.empty()))
  {
    fprintf(stderr, "Missing dimacs/population/distance or database.\n");
    exit(1);
  }

  if(!database.empty() && level.empty())
  {
    fprintf(stderr, "Missing level for database.\n");
    exit(1);
  }

  if(!seen_state && strlen(state) < 2)
  {
    fprintf(stderr, "Missing state.\n");
    exit(1);
  }

  if(!database.empty())
  {
    //FIXME use filesystem when switched to C++17 or greater
#ifdef _WIN32
    char sep = '\\';
#else
    char sep = '/';
#endif
    if(level == "tracts") level = "census_tracts";
    else if(level != "counties")
    {
      fprintf(stderr, "Unknown level %s.\n", level.c_str());
      exit(1);
    }
    rp.dimacs_file     = database + sep + rp.state + sep + level + sep + "graph" + sep + rp.state + ".dimacs";
    rp.population_file = database + sep + rp.state + sep + level + sep + "graph" + sep + rp.state + ".population";
    rp.distance_file   = database + sep + rp.state + sep + level + sep + "graph" + sep + rp.state + "_distances.csv";
  }


  // debug
  cout << "dimacs_file     = " << rp.dimacs_file << endl;
  cout << "distance_file   = " << rp.distance_file << endl;
  cout << "population_file = " << rp.population_file << endl;
  cout << "state[2]        = " << rp.state << endl;
  cout << "L               = " << rp.L << endl;
  cout << "U               = " << rp.U << endl;
  cout << "k               = " << rp.k << endl;
  cout << "model           = " << rp.model << endl;
  cout << "ralg_hot_start  = " << rp.ralg_hot_start << endl;
  cout << "output          = " << rp.output << endl;

  return rp;
}
