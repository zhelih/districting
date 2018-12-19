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
      sscanf(buf, "p edge %d", &nr_nodes);
//      sscanf(buf, "p edge %d %d", &nr_nodes, &nr_edges);
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

  for(int i = 0; i < nr_nodes; ++i)
    sort(g->nb(i).begin(), g->nb(i).end());

  fclose(f);
  return g;
}

graph::graph(uint n)
{
  nr_nodes = n;
  nb_.resize(n);
}

graph::~graph()
{
}

void graph::add_edge(uint i, uint j)
{
  if(is_edge(i,j))
    return;
  nb_[i].push_back(j);
  nb_[j].push_back(i);
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

bool graph::is_edge(uint i, uint j)
{
  for(uint k : nb_[i])
    if(k == j)
      return true;
  for(uint k : nb_[j])
    if(k == i)
      return true;
  return false;
}
