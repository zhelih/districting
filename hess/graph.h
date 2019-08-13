#ifndef _GRAPH_H
#define _GRAPH_H

typedef unsigned int uint;

#include <vector>
#include <stack>

using namespace std;

class graph
{
private:
    std::vector<std::vector<int> > nb_;
    bool is_edge(uint i, uint j);
    int k;
public:
    uint nr_nodes;
    graph(uint n);
    virtual ~graph();
    void add_edge(uint i, uint j);
    void remove_edge(uint i, uint j);
    std::vector<int>& nb(uint i) { return nb_[i]; }
    bool is_connected(); // TODO const;

    // works as far as no pointers are members
    graph* duplicate() const { return new graph(*this); }
    int get_k() const;
    void set_k(int k_) { k = k_; }
    void connect(const std::vector<std::vector<int>>& dist); // make the graph connected
};

graph* from_dimacs(const char* fname); // don't forget to delete

//FIXME
#ifndef INFINITY
#define INFINITY 1e20
#endif



#endif
