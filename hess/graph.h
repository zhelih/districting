#ifndef _GRAPH_H
#define _GRAPH_H

typedef unsigned int uint;

#include <vector>

class graph
{
  private:
  std::vector<std::vector<int> > nb_;
  bool is_edge(uint i, uint j) ;
  public:
  uint nr_nodes;
  graph(uint n);
  ~graph();
  void add_edge(uint i, uint j);
  std::vector<int>& nb(uint i) { return nb_[i]; }
  bool is_connected(); // TODO const;
};

graph* from_dimacs(const char* fname); // don't forget to delete

#endif
