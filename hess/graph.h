#ifndef _GRAPH_H
#define _GRAPH_H

typedef unsigned int uint;

#include <vector>
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
    vector<int> preprocess(vector<int>& new_population, vector<bool>& deleted, int L, int U);
    vector<int> findUnderPopulatedLeaves(vector<int> new_population, vector<bool> deleted, int L);
    void clean(vector<int>& new_population, vector<bool>& deleted, int L, int U, int& numOfEdgeDel, int& numOfNodeMerge);
    void edgeClean(vector<int>& population, int U);
    //bool graph::deleteEdge(int i, int j);
    // works as far as no pointers are members
    graph* duplicate() const { return new graph(*this); }
    int get_k() const { return k; }
    void set_k(int k_) { k = k_; }
};

graph* from_dimacs(const char* fname); // don't forget to delete

#endif
