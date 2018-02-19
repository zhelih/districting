#include "graph.h"
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

//#define FROM_1

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

graph* from_dimacs(const char* fname)
{
  FILE* f = fopen(fname, "r");
  if(!f)
  {
    fprintf(stderr, "Cannot open file %s\n", fname);
    return 0;
  }

  char buf[1024]; int nr_nodes = 0, nr_edges = 0;

  // search for "n nodes edges"
  while(NULL != fgets(buf, sizeof(buf), f))
  {
    if(buf[0] == 'p')
    {
      sscanf(buf, "p edge %d %d", &nr_nodes, &nr_edges);
      break;
    }
  }

  if(nr_nodes == 0 && nr_edges == 0)
  {
    fclose(f);
    fprintf(stderr, "Cannot found dimacs metadata in %s\n", fname);
    return 0;
  }

//  fprintf(stderr, "Found metadata in %s : (%d, %d)\n", fname, nr_nodes, nr_edges);

  graph* g = new graph(nr_nodes);

  rewind(f);

  //read the edges
  while(NULL != fgets(buf, sizeof(buf), f))
  {
    if(buf[0] == 'e')
    {
      int from, to; double weight;
      sscanf(buf, "e %d %d %lf", &from, &to, &weight);
#ifdef FROM_1
      from--; to--;
#endif
      g->add_edge(from, to);//, weight);
    }
  }

  printf("graph: %d nodes, %d edges\n", nr_nodes, nr_edges);

  fclose(f);
  return g;
}

graph::graph(uint n)
{
  // compute bit mask
  mask[0] = 1;
  for(uint i = 1; i < CHUNK_SIZE; ++i)
    mask[i] = mask[i-1] << 1;

 // create adj matrix
  //adj = new uint*[n];
  adj.resize(n);
  for(uint i = 0; i < n; ++i)
  {
    //adj[i] = new uint[n];
    adj[i].resize(n);
  }
  for(uint i = 0; i < n; ++i)
    for(uint j = 0; j < n; ++j)
      adj[i][j] = 0;
  nr_nodes = n;
  nb_.resize(n);
}

graph::~graph()
{
/*  for(uint i = 0; i < nr_nodes; ++i)
    delete [] adj[i];
  delete [] adj;*/
}

void graph::add_edge(uint i, uint j)
{
  if(!adj[i][j]) nb_[i].push_back(j); // prevent multiple edges
  if(!adj[j][i]) nb_[j].push_back(i);
  adj[i][j] = true;
  adj[j][i] = true;
}

bool graph::is_connected()
{
  vector<int> s;
  s.push_back(0);
  int *visited = new int[nr_nodes];
  for(uint i = 0; i < nr_nodes; ++i)
    visited[i] = 0;
  while(!s.empty())
  {
    int cur = s.back(); s.pop_back();
    visited[cur] = 1;
    auto nbs = nb(cur);
    for(auto it = nbs.begin(); it != nbs.end(); ++it)
    {
      if(!visited[*it])
      {
        s.push_back(*it);
      }
    }
  }
  bool res = true;
  for(unsigned int i = 0; i < nr_nodes; ++i)
    if(!visited[i])
    {
      res = false;
      break;
    }
  delete [] visited;
  return res;
}

/*void graph::check_solution()
{
  for(i = 0; i < n; ++i)
    if(x[ii] == 1)
    {
      for(j=i+1; j < n; ++j)
      {
        if(x[jj] == 1)
        {
          if(! i and j connected via xia and xjb)
          {
            compute min i,j-separator out of xia and xjb
            add lazy
          }
        }
      }
    }
}
*/
/*inline bool graph::is_edge(uint i, uint j)
{
  return adj[i][j/CHUNK_SIZE]&mask[j%CHUNK_SIZE];
}*/
