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
    
    void edgeClean(const vector<int>& population, int U);
    void edgeCleanNeighbor(vector<int>& new_population, int s, int U);
    vector<vector<int>> FindBiconnectedComponents(vector<int> &AV, vector<bool> &deletedNodes);
    std::vector<std::vector<int>> FindConnectedComponents(vector<bool> &deletedNodes);
    void Bico_Sub(int v, int u, int &i, vector<int> &number, vector<int> &lowopt, stack<int> &le, stack<int> &re, vector< vector<int>> &BC, vector<bool> &deletedNodes);
    vector<int> ShortestPathsNodeWeighted(int origin, const vector<int>& population, vector<bool> &deletedNodes, int total_pop);
    //bool graph::deleteEdge(int i, int j);
    // works as far as no pointers are members
    graph* duplicate() const { return new graph(*this); }
    int get_k() const { return k; }
    void set_k(int k_) { k = k_; }
    void connect(const std::vector<std::vector<int>>& dist); // make the graph connected    
};

graph* from_dimacs(const char* fname); // don't forget to delete

#endif