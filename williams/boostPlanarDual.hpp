#ifndef __CREATE_DUAL_GRAPH_HPP__
#define __CREATE_DUAL_GRAPH_HPP__

#include <vector>
#include <fstream>
#include <iostream>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/planar_face_traversal.hpp>

std::ofstream myfile ("primalGraph.txt");
long s = 0;

namespace boost

{

  template <typename InputGraph,
            typename OutputGraph,
            typename EdgeIndexMap>
  struct dual_graph_visitor : public planar_face_traversal_visitor
  {
    //typedef typename myfile;
    typedef typename graph_traits<OutputGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<InputGraph>::edge_descriptor edge_t;
    typedef typename std::vector<vertex_t> vertex_vector_t;
    typedef iterator_property_map
      < typename vertex_vector_t::iterator, EdgeIndexMap >
        edge_to_face_map_t;

    dual_graph_visitor(InputGraph& arg_g,
                       OutputGraph& arg_dual_g,
                       EdgeIndexMap arg_em
                       ) :
      g(arg_g),
      dual_g(arg_dual_g),
      em(arg_em),
      edge_to_face_vector(num_edges(g),
                          graph_traits<OutputGraph>::null_vertex()),
      edge_to_face(edge_to_face_vector.begin(), em)
      {}

    void begin_face()
    {
      current_face = add_vertex(dual_g);
    }

    template <typename Edge>
    void next_edge(Edge e)
    {
      if (s==0)
      {
          myfile << "vertices: " << boost::num_vertices(g) << std::endl;
          myfile << "edges: " << boost::num_edges(g) << std::endl;
      }
      s++;
      vertex_t existing_face = edge_to_face[e];
      if (existing_face == graph_traits<OutputGraph>::null_vertex())
      {
        edge_to_face[e] = current_face;
      }
      else
      {
        add_edge(existing_face, current_face, dual_g);
        myfile << e << std::endl;
      }
    }

    InputGraph& g;
    OutputGraph& dual_g;
    EdgeIndexMap em;
    vertex_t current_face;
    vertex_vector_t edge_to_face_vector;
    edge_to_face_map_t edge_to_face;
  };

        template <typename InputGraph,
            typename OutputGraph,
            typename PlanarEmbedding,
            typename EdgeIndexMap
           >
  void create_dual_graph(InputGraph& g,
                         OutputGraph& dual_g,
                         PlanarEmbedding embedding,
                         EdgeIndexMap em)
  {
    dual_graph_visitor<InputGraph, OutputGraph, EdgeIndexMap>
    visitor(g, dual_g, em);
    planar_face_traversal(g, embedding, visitor, em);
  }


  template <typename InputGraph,
            typename OutputGraph,
            typename PlanarEmbedding
           >
  void create_dual_graph(InputGraph& g,
                         OutputGraph& dual_g,
                         PlanarEmbedding embedding
                         )
  {
    create_dual_graph(g, dual_g, embedding, get(edge_index,g));
  }

}

#endif //__CREATE_DUAL_GRAPH_HPP__
