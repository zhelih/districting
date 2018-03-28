#ifndef _GRAPH_H
#define _GRAPH_H

typedef unsigned int uint;

#include <vector>

#define CHUNK_SIZE (8*(sizeof(int))) // chunk size
class graph
{
  private:
  std::vector<std::vector<int> >  adj;
  std::vector<std::vector<int> > nb_;
  uint mask[CHUNK_SIZE];
  public:
  uint nr_nodes;
  graph(uint n);
  ~graph();
  void add_edge(uint i, uint j);
  inline bool is_edge(uint i, uint j) {return adj[i][j];}
  std::vector<int>& nb(uint i) { return nb_[i]; }
  bool is_connected(); // TODO const;
  void floyd_warshall(std::vector<std::vector<int>>& d);
//  inline uint weight(uint i) { return 1; }
};

graph* from_dimacs(const char* fname); // don't forget to delete
graph* from_grid(int n);

#endif
