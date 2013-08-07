
#include <point-process-core/marked_grid.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <math-core/geom.hpp>
#include <planner-core/edmonds_optimum_branching.hpp>
#include <iostream>
#include <iterator>


using namespace boost;
using namespace point_process_core;
using namespace math_core;


int main()
{

  // the marked grid of visited nodes
  nd_aabox_t window;
  window.start = point( 0 );
  window.end = point( 10 );
  double res = 1.0;
  marked_grid_t<bool> visisted_grid( window, res );

  // the current position
  nd_point_t current_position = point( 5.5 );

  // types of the graph
  typedef adjacency_list < listS, vecS, directedS,
			   no_property, property < edge_weight_t, float > > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
  typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
  typedef std::pair<int,int> Edge;
  
  // calculate the current cell
  marked_grid_cell_t current_cell = visisted_grid.cell( current_position );
  
  // create the edge set
  std::vector<Edge> edge_set;
  std::vector<float> weight_set;
  long num_nodes = 0;
  int current_cell_vertex_number = -1;
  std::vector<marked_grid_cell_t> all_cells 
    = visisted_grid.all_cells();
  for( size_t i = 0; i < all_cells.size(); ++i ) {
    
    // skip already visited cells
    // but hte currrent position is special since it is our source
    if( visisted_grid( all_cells[i] ) == true &&
	all_cells[i] != current_cell )
      continue;

    // make curernt cell number 
    if( all_cells[i] == current_cell ) {
      current_cell_vertex_number = i;
    }
    
    // increase the number of nodes count
    ++num_nodes;
    
    // ok, add and edge from this cell to all other non-visited cells
    // with the ditance of the centers of the cells as the edge weight
    for( size_t k = 0; k < all_cells.size(); ++k ) {
      if( visisted_grid( all_cells[k] ) == true ||
	  i == k)
	break;
      
      // Calculate distance from cells
      nd_aabox_t start_region = visisted_grid.region( all_cells[i] );
      nd_aabox_t end_region = visisted_grid.region( all_cells[k] );
      nd_point_t start_point = start_region.start + ( start_region.end - start_region.start ) * 0.5;
      nd_point_t end_point = end_region.start + ( end_region.end - end_region.start ) * 0.5;
      float dist = distance( start_point, end_point );
      
      // add the edge (fowards and backwards)
      edge_set.push_back( Edge( i, k ) );
      weight_set.push_back( dist );
      edge_set.push_back( Edge( k, i ) );
      weight_set.push_back( dist );
      
    }
  }

  // add a single edge from a new "root" node to the current cell
  int root_vertex_number = current_cell_vertex_number;
  // int root_vertex_number = num_nodes;
  // ++num_nodes;
  // edge_set.push_back( Edge( root_vertex_number, current_cell_vertex_number ) );
  // weight_set.push_back( 1.0 );
  // edge_set.push_back( Edge( current_cell_vertex_number, root_vertex_number ) );
  // weight_set.push_back( 1.0 );

  
  // Create a graph from the edge and weight set
  graph_t graph( edge_set.begin(), edge_set.end(), weight_set.begin(), num_nodes );
  
  // get the vertex index map and weight map
  property_map<graph_t,edge_weight_t>::type weightmap = get(edge_weight,graph);
  property_map<graph_t,vertex_index_t>::type indexmap = get(vertex_index, graph);

  // get the wanted roots fo the branchings
  vertex_descriptor s = vertex( root_vertex_number, graph );
  std::vector<vertex_descriptor> roots;
  roots.push_back( s );
  
  // ok, get a storage for the branchings
  std::vector<edge_descriptor> branchings;
  
  // find hte branchings of the graph
  edmonds_optimum_branching<false,true,true>
    ( graph,
      indexmap,
      weightmap,
      roots.begin(),
      roots.end(),
      std::back_inserter(branchings) );
  
  // print out the branchings
  for( size_t i = 0; i < branchings.size(); ++i ) {
    std::cout << "[" << i << "]  " << branchings[i] << std::endl;
  }
  
  return 0;
}
