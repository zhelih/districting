#include <iostream>
#include <vector>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <C:\Users\hvalidi\Desktop\Redistricting\planar_dual.hpp>
using namespace boost;
using namespace std;

template <typename Graph>

void print_graph(Graph& g)
{
  ofstream myfile ("dualGraph.txt");

  typename graph_traits<Graph>::vertex_iterator vi, vi_end;

  myfile << "vertices: " << boost::num_vertices(g) << std::endl;

  typename graph_traits<Graph>::edge_iterator ei, ei_end;

  myfile << "edges: " << boost::num_edges(g) << std::endl;

  for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
  {
    myfile << *ei << std::endl;

  }

}

int main(int argc, char** argv)
{
   typedef adjacency_list
    < vecS,
     vecS,
     undirectedS,
     property<vertex_index_t, int>,
     property<edge_index_t, int>
    >
    graph;

    // Create the graph - two triangles that share an edge
    int a,b,m,n;
    string c, temp;
    ifstream input("TX.txt");
    input >> temp >> temp >> n >> m;
    graph g(n);

    for(int i=0; i < m; i++)
    {
        input >> temp >> a >> b;
        add_edge(a,b,g);
    }

    vector<int> listOfEdges (m);

    // Create an empty graph to hold the dual
    graph dual_g;

    property_map<graph, edge_index_t>::type e_index = get(edge_index, g);

    graph_traits<graph>::edges_size_type edge_count = 0;

    graph_traits<graph>::edge_iterator ei, ei_end;

    for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    {
      put(e_index, *ei, edge_count++);
    }

    // Compute the planar embedding
    typedef std::vector< graph_traits<graph>::edge_descriptor > vec_t;
    std::vector<vec_t> embedding(num_vertices(g));
    bool planar = boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                               boyer_myrvold_params::embedding = &embedding[0]
                               );

    cerr<<"Is graph planar? "<<planar<<endl;

    create_dual_graph(g, dual_g, &embedding[0]);

    print_graph(dual_g);

  return 0;
}
