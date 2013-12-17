
#include "shortest_path_next_planner.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include <math-core/utils.hpp>
#include "edmonds_optimum_branching.hpp"
#include <point-process-core/context.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <plot-server/api/plot.hpp>
#include <boost/property_tree/ptree.hpp>


using namespace point_process_core;
using namespace boost::property_tree;
using namespace plot_server::api;

namespace planner_core {

  //=======================================================================

  std::vector<marked_grid_cell_t>
  canonical_point_bin( const std::vector<nd_point_t>& points,
		       const marked_grid_t<bool>& grid )
  {
    for( size_t i = 0; i < points.size(); ++i ) {
      std::cout << points[i] << " ";
    }
    std::cout << std::endl;
    std::vector<marked_grid_cell_t> bins;
    std::vector<marked_grid_cell_t> cells = grid.all_cells();
    for( size_t i = 0; i < cells.size();++i ) {
      //bool found = false;
      for( size_t k = 0; k < points.size(); ++k ) {
	nd_aabox_t region = grid.region( cells[i] );
	if( is_inside( points[k], region ) ) {
	  //found = true;
	  bins.push_back( cells[i] );
	  break;
	}
      }
    }

    std::sort( bins.begin(),
	       bins.end() );

    return bins;
  }

  //=======================================================================

  void spnp_console_print_point_set( const std::vector<nd_point_t>& points,
				     const marked_grid_t<bool>& grid )
  {
    std::vector<marked_grid_cell_t> full_cells;
    for( size_t i = 0; i < points.size(); ++i ) {
      marked_grid_cell_t cell = grid.cell( points[i] );
      if( std::find( full_cells.begin(),
		     full_cells.end(),
		     cell ) == full_cells.end() ) {
	full_cells.push_back( cell );
      }
    }
    
    // Ok, now get all possible cells
    std::vector<marked_grid_cell_t> all_cells = grid.all_cells_ordered();
    
    // Print out the cells, in a grid format
    std::vector<std::string> barrier_characters;
    barrier_characters.push_back( "" );
    barrier_characters.push_back( "|" );
    barrier_characters.push_back( "*" );
    int barrier_characters_len = barrier_characters.size();
    std::vector<std::string> end_characters;
    end_characters.push_back( "\n" );
    int end_characters_len = end_characters.size();
    char join_barrier_character = '+';
    std::string join_end_string = "\n";
    int barrier_spacing = 5;
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      marked_grid_cell_t cell = all_cells[i];
      bool marked = (std::find( full_cells.begin(),
				full_cells.end(),
				cell ) != full_cells.end());
      
      // write out empoty of filled cell
      if( marked ) {
	std::cout << "v";
      } else {
	if( grid( cell ) ) {
	  std::cout << ".";
	} else {
	  std::cout << "-";
	}
      }
      
      // write out a barrier if needed,
      // first track how many barriers are reached (count) and 
      // the last index for a barrier
      int barrier_count = 0;
      int barrier_idx = 0;
      int end_count = 0;
      int end_idx = 0;
      nd_point_t start_point = centroid( grid.region( all_cells[0] ) );
      nd_point_t cell_point = centroid( grid.region( cell ) );
      for( int k = 0; k < start_point.n; ++k ) {
	double diff = abs(start_point.coordinate[k] - cell_point.coordinate[k]);
	double cell_size = grid.cell_sizes()[k];
	int grid_diff = (int)symmetric_round( diff / cell_size );
	
	if( (grid_diff + 1) % barrier_spacing == 0 && grid_diff > 0 ) {
	  barrier_count++;
	  barrier_idx = k;
	}
	nd_vector_t dir = zero_vector( start_point.n );
	dir.component[k] = 0.5 * cell_size;
	nd_vector_t all_dirs = 0.5 * (point( grid.cell_sizes() ) - zero_point( start_point.n ) );
	
	if( k > 0 &&
	    is_inside( cell_point, grid.window() ) == false &&
	    is_inside( cell_point + (-1.0 * dir), grid.window() ) ) {
	  ++end_count;
	  end_idx = k;
	} 
      }
      
      // Ok, use join barrier if more that one reached, else
      // just use the indexed character for a single barrier
      if( barrier_count > 1 ) {
	std::cout << join_barrier_character;
      } else if( barrier_count == 1 ) {
	std::cout << barrier_characters[ barrier_idx % barrier_characters_len ];
      }
      
      // Write out end character if needed
      if( end_count > 1 ) {
	std::cout << join_end_string;
      } else if( end_count == 1 ) {
	std::cout << end_characters[ end_idx % end_characters_len ];
      }
    }
    std::cout << std::endl;
    std::cout.flush();
  }
  
  //=======================================================================

  void spnp_console_print_next_obs_dist( const marked_grid_t<double>& grid )
  {
    for( double xi = 0.5; xi < 11; xi += 1 ) {
      boost::optional<double> mark = grid( point(xi) );
      std::cout << "|";
      std::cout << std::setw(4) << std::setprecision(4);
      if( mark ) {
	std::cout << *mark;
      } else {
	std::cout << "0";
      }
    }
    std::cout << "|";
    std::cout << std::endl;
    std::cout.flush();
  }

  //=======================================================================

  // boost::optional<marked_grid_cell_t>
  // shortest_path_next_planner
  // ::next_observation_in_shortest_path
  // ( const std::vector<nd_point_t>& points,
  //   const marked_grid_t<bool> visisted_grid_param,
  //   const nd_point_t& current_position ) const
  // {

  //   // This is a hacked version which assumes we are going from left to right
  //   // and so the next shortest is the first non-visited cell from
  //   // left to right
    
  //   // append to visited cells any cells without points in them
  //   marked_grid_t<bool> visited_grid = visisted_grid_param.copy_structure<bool>();
  //   std::vector<marked_grid_cell_t> cells = visisted_grid_param.all_marked_cells();
  //   for( size_t i = 0; i < cells.size(); ++i ) {
  //     visited_grid.set( cells[i], true );
  //   }
  //   cells = visited_grid.all_cells();
  //   for( size_t k = 0; k < cells.size(); ++k ) {
  //     nd_aabox_t region = visited_grid.region( cells[k] );
  //     bool found = false;
  //     for( size_t i = 0; i < points.size(); ++i ) {
  // 	nd_point_t p = points[i];
  // 	if( is_inside( p, region ) ) {
  // 	  found = true;
  // 	  break;
  // 	}
  //     }
  //     if( found == false ) {
  // 	visited_grid.set( cells[k], true );
  //     }
  //   }
   
  //   // Ok, now get all unvisited cells and take their starting points
  //   // diostance to current location
  //   std::vector<double> distances;
  //   std::vector<marked_grid_cell_t> cell_for_dist;
  //   cells = visited_grid.all_cells();
  //   for( size_t i = 0; i < cells.size(); ++i ) {
  //     if( visited_grid( cells[i] ) ) {
  // 	continue;
  //     }
  //     nd_aabox_t region = visited_grid.region( cells[i] );
  //     distances.push_back( distance_sq( region.start, current_position ) );
  //     cell_for_dist.push_back( cells[i] );
  //   }

  //   // return nothing if no possible cells
  //   if( distances.empty() ) {
  //     return boost::optional<marked_grid_cell_t>();
  //   }

  //   // find the minimum starting point, this is the
  //   // next best action
  //   std::vector<double>::iterator min_iter
  //     = std::min_element( distances.begin(),
  // 			  distances.end() );
  //   size_t min_idx = std::distance( distances.begin(), min_iter );
    
  //   // return hte cell for the min distance
  //   return boost::optional<marked_grid_cell_t>( cell_for_dist[ min_idx ] );
  // }

  //=======================================================================

  boost::optional<marked_grid_cell_t>
  shortest_path_next_planner
  ::next_observation_in_shortest_path
  ( const std::vector<nd_point_t>& points,
    const marked_grid_t<bool> visisted_grid_param,
    const nd_point_t& current_position ) const
  {
    
    // we are going to make a graph of all of the possible
    // observation regions/cells left and have edges with the
    // path distance to them, then find the shortest path
    // using dijkstras
    
    // make use of boost graph without needing the boost;:
    using namespace boost;
    
    // types of the graph
    typedef adjacency_list < listS, vecS, directedS,
  			     no_property, property < edge_weight_t, float > > graph_t;
    typedef graph_traits < graph_t >::vertex_descriptor vertex_descriptor;
    typedef graph_traits < graph_t >::edge_descriptor edge_descriptor;
    typedef std::pair<int,int> Edge;
    
    // calculate the current cell
    marked_grid_cell_t current_cell = visisted_grid_param.cell( current_position );

    // append to visited cells any cells without points in them
    marked_grid_t<bool> visited_grid = visisted_grid_param.copy_structure<bool>();
    std::vector<marked_grid_cell_t> cells = visisted_grid_param.all_marked_cells();
    for( size_t i = 0; i < cells.size(); ++i ) {
      visited_grid.set( cells[i], true );
    }
    cells = visited_grid.all_cells();
    for( size_t k = 0; k < cells.size(); ++k ) {
      nd_aabox_t region = visited_grid.region( cells[k] );
      bool found = false;
      for( size_t i = 0; i < points.size(); ++i ) {
  	nd_point_t p = points[i];
  	if( is_inside( p, region ) ) {
  	  found = true;
  	  break;
  	}
      }
      if( found == false ) {
  	visited_grid.set( cells[k], true );
      }
    }
    
    // create the edge set
    std::vector<Edge> edge_set;
    std::vector<float> weight_set;
    std::vector<marked_grid_cell_t> vertex_index_to_cell;
    std::map<marked_grid_cell_t,int> cell_to_vertex_index;
    long num_nodes = 0;
    int current_cell_vertex_number = -1;
    std::vector<marked_grid_cell_t> all_cells 
      = visited_grid.all_cells();
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      
      // skip already visited cells
      // but hte currrent position is special since it is our source
      if( visited_grid( all_cells[i] ) == true &&
  	  all_cells[i] != current_cell )
  	continue;
      
      
      // make curernt cell number 
      if( all_cells[i] == current_cell ) {
  	current_cell_vertex_number = i;
      }
      
      // increase the number of nodes count
      if( cell_to_vertex_index.find( all_cells[i] ) ==
  	  cell_to_vertex_index.end() ) {
  	vertex_index_to_cell.push_back( all_cells[i] );
  	cell_to_vertex_index[ all_cells[i] ] = num_nodes;
  	++num_nodes;
      }	
      
      // ok, add and edge from this cell to all other non-visited cells
      // with the ditance of the centers of the cells as the edge weight
      for( size_t k = 0; k < all_cells.size(); ++k ) {
  	if( ( visited_grid( all_cells[k] ) == true &&
  	      all_cells[k] != current_cell ) ||
  	    i == k)
  	  continue;

  	// update cell-to-index 
  	if( cell_to_vertex_index.find( all_cells[k] ) ==
  	    cell_to_vertex_index.end() ) {
  	  vertex_index_to_cell.push_back( all_cells[k] );
  	  cell_to_vertex_index[ all_cells[k] ] = num_nodes;
  	  ++num_nodes;
  	}	

  	// Calculate distance from cells
  	nd_aabox_t start_region = visited_grid.region( all_cells[i] );
  	nd_aabox_t end_region = visited_grid.region( all_cells[k] );
  	nd_point_t start_point = start_region.start + ( start_region.end - start_region.start ) * 0.5;
  	nd_point_t end_point = end_region.start + ( end_region.end - end_region.start ) * 0.5;
  	float dist = distance( start_point, end_point );
	
  	// add the edge (fowards and backwards)
  	edge_set.push_back( Edge( cell_to_vertex_index[all_cells[i]], 
  				  cell_to_vertex_index[all_cells[k]] ) );
  	weight_set.push_back( dist );
  	//edge_set.push_back( Edge( k, i ) );
  	//weight_set.push_back( dist );
	
      }
    }
    
    // add a single edge from a new "root" node to the current cell
    current_cell_vertex_number = cell_to_vertex_index[ current_cell ];
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
    
    // so this branching is 2*optimal at worst, so 
    // we will use it!
    
    // resturn the cell for the first edge in the branching
    // or nothing if there were no branchings found
    if( branchings.empty() ) {

      // debug 
      //std::cout << "      edmonds found no roots" << std::endl;
  
      return boost::optional<marked_grid_cell_t>();
    } else {
      return boost::optional<marked_grid_cell_t>(vertex_index_to_cell[ target(branchings[0], graph) ]);
    }
  }

  //=======================================================================

  void
  shortest_path_next_planner::
  estimate_shortest_path_next_observation_distribution
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const marked_grid_t<bool>& visited_grid,
    const unsigned long num_samples_of_point_sets,
    const nd_point_t& current_position,
    marked_grid_t<double>& next_observation_distribution) const
  {

    scoped_context_switch context( chain_context( context_t( "next-obs-dist" ) ) );

    // We will sample from the process some points sets.
    // Then, foir each we will estiamte hte next location in the shortest
    // path through the points.
    // We will then increment the count of that location
    // At the end we will divide counts by total samples to get
    // an estiamte of the distribution oif best nex observation
    
    // take samples from the process
    for( size_t sample_i = 0; sample_i < num_samples_of_point_sets; ++sample_i ) {
      
      // sample a point set
      std::vector<nd_point_t> sample_point_set 
	= process->sample_and_step();

      // skip updating distribution if no points in sample
      if( sample_point_set.empty() ) {

	// debug 
	std::cout << " #0 ";
	std::cout.flush();

	continue;
      }

      // print this sample to the console
      spnp_console_print_point_set( sample_point_set, visited_grid );
      
      
      // Calcualte the next observation in the shortest path thorugh
      // points
      boost::optional<marked_grid_cell_t> next_obs 
	= next_observation_in_shortest_path( sample_point_set,
					     visited_grid,
					     current_position );
      
      // skip updating distribution if no next observation found
      if( !next_obs ) {
	std::cout << std::endl;
	continue;
      }
      
      // Ok, increment the mark at this cell
      boost::optional<double> m = next_observation_distribution( *next_obs );
      if( m ) {
	next_observation_distribution.set( *next_obs, *m + 1.0 );
      } else {
	next_observation_distribution.set( *next_obs, 1.0 );
      }

      // print out the chosen next obs and hte ddistribution counts
      std::cout << std::setw(5) << *next_obs << "  ";
      spnp_console_print_next_obs_dist( next_observation_distribution );
    }

    // // Ok, now divide all counts by number of samples
    // std::vector<marked_grid_cell_t> cells = next_observation_distribution.all_marked_cells();
    // for( size_t i = 0; i < cells.size(); ++i ) {
    //   next_observation_distribution.set( cells[i],
    // 					 *next_observation_distribution( cells[i] ) / (double)num_samples_of_point_sets );
    // }

  }


  //=======================================================================

  size_t nonvisited_grid_count( const std::vector<nd_point_t>& sample, 
				const marked_grid_t<bool>& visited_grid )
  {
    std::vector<marked_grid_cell_t> new_cells;
    for( auto p : sample ) {
      marked_grid_cell_t cell = visited_grid.cell( p );
      if( !visited_grid( cell ) ) {
	if( std::find( new_cells.begin(), new_cells.end(), cell ) 
	    == new_cells.end() ) {
	  new_cells.push_back( cell );
	}
      }
    }
    return new_cells.size();
  }
  
  //=======================================================================

  std::vector<marked_grid_cell_t>
  estimate_most_likely_world
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const marked_grid_t<bool>& visited_grid,
    const unsigned long num_samples_of_point_sets,
    const unsigned long skip_between_samples_used = 10 )
  {
    
    // a list of bins and their counts
    std::vector<double> bin_counts;
    std::vector< std::vector<marked_grid_cell_t> > bins;


    // take samples from the process
    for( size_t sample_i = 0; sample_i < num_samples_of_point_sets; ++sample_i ) {
      
      // sample a point set
      std::vector<nd_point_t> sample_point_set 
	= process->sample_and_step();
      

      // keep sampling until we get a point set that has at least some
      // possible new information
      // (so, at least one *new* non-visited grid must have a point
      size_t max_find_samples = 100;
      size_t find_samples_count = 0;
      while( nonvisited_grid_count( sample_point_set, visited_grid ) < 1
	     && find_samples_count < max_find_samples ) {
	sample_point_set = process->sample_and_step();
	++find_samples_count;
      }

      if( nonvisited_grid_count( sample_point_set, visited_grid ) == 0 ) {
	std::cout << "  !! no new cells in " << max_find_samples << " samples..." << std::endl;
      }

      std::vector<marked_grid_cell_t> sample_bin = canonical_point_bin( sample_point_set, visited_grid );

      // skip the given samples
      for( unsigned long skip_i = 0; skip_i + 1 < skip_between_samples_used; ++skip_i ) {
	process->single_mcmc_step();
      }

      // skip if no points in sample
      if( sample_point_set.empty() ) {

	// debug 
	std::cout << " #0 ";
	std::cout.flush();

	continue;
      }

      // print this sample to the console
      spnp_console_print_point_set( sample_point_set, visited_grid );
      std::cout << std::endl;

      // print the canonical representation
      if( false ) {
	for( size_t i = 0; i < sample_bin.size(); ++i ) {
	  std::cout << sample_bin[i];
	}
	std::cout << std::endl;
      }

      // ok, increment count of bin
      std::vector< std::vector<marked_grid_cell_t> >::iterator fiter
	= std::find( bins.begin(),
		     bins.end(),
		     sample_bin );
      if( fiter == bins.end() ) {
	bins.push_back( sample_bin );
	bin_counts.push_back( 1 );
      } else {
	size_t fiter_idx = std::distance( bins.begin(), fiter );
	bin_counts[ fiter_idx ] += 1;
      }
    }


    // print out the table of worlds
    std::cout << "Worlds: " << std::endl;
    for( size_t i = 0; i < bins.size(); ++i ) {
      std::cout << "  [" << std::setw(4) << bin_counts[i] << "] ";
      std::vector<marked_grid_cell_t> world = bins[i];
      for( size_t k = 0; k < world.size(); ++k ) {
	std::cout << "(";
	for( size_t r = 0; r < world[k].coordinate.size(); ++r ) {
	  std::cout << world[k].coordinate[r];
	  if( r + 1 < world[k].coordinate.size() ) {
	    std::cout << ",";
	  }
	}
	std::cout << ")";
      }
      std::cout << std::endl;
    }

    // Ok, now find the maximum count bin and return it
    std::vector< double >::iterator max_iter 
      = std::max_element( bin_counts.begin(),
			  bin_counts.end() );
    size_t max_idx = std::distance( bin_counts.begin(), max_iter );

    // print the best world
    std::cout << "Most Likely world: ";
    std::vector<marked_grid_cell_t> world = bins[max_idx];
    for( size_t k = 0; k < world.size(); ++k ) {
      std::cout << world[k];
    }
    std::cout << std::endl;

    return bins[ max_idx ];
  }

  // //=======================================================================

  // marked_grid_cell_t 
  // shortest_path_next_planner::
  // choose_next_observation_cell()
  // {

  //   scoped_context_switch context( chain_context( context_t( "greedy-shortest-path-choose-next-cell" ) ) );
    
  //   // For all non-visited cells,
  //   // compute the distribution over the expected
  //   // next action for shortest path to all points
  //   // given we take the observation
  //   marked_grid_t<double> next_observation_distribution = _visited_grid.copy_structure<double>();
  //   std::vector<marked_grid_cell_t> all_cells 
  //     = _visited_grid.all_cells();

  //   // estimate the next observation distribution given an observation
  //   estimate_shortest_path_next_observation_distribution
  //     ( _point_process,
  // 	_visited_grid,
  // 	_sampler_planner_params.num_samples_of_point_sets,
  // 	_current_position,
  // 	next_observation_distribution );

  //   // print out the coutns
  //   spnp_console_print_next_obs_dist( next_observation_distribution );

  //   // normalize the cells to be probabilities
  //   all_cells = next_observation_distribution.all_marked_cells();
  //   double sum = 0;
  //   for( size_t i = 0; i < all_cells.size(); ++i ) {
  //     sum += *next_observation_distribution(all_cells[i]);
  //   }
  //   for( size_t i = 0; i < all_cells.size(); ++i ) {
  //     next_observation_distribution.set( all_cells[i], *next_observation_distribution(all_cells[i]) / sum );
  //   }

  //   // Ok, now get the set of actions which have a high enough
  //   // probability
  //   std::vector<marked_grid_cell_t> cells = next_observation_distribution.all_marked_cells();
  //   std::vector<marked_grid_cell_t> potential_next_obs;
  //   for( size_t i = 0; i < cells.size(); ++i ) {
  //     marked_grid_cell_t cell = cells[i];
  //     if( *next_observation_distribution(cell) >= _probability_next_action_is_in_shortest_path_threshold ) {
  // 	potential_next_obs.push_back( cell );
  //     }
  //   }

  //   // if we have no potentials, add all the marked cells
  //   if( potential_next_obs.empty() ) {

  //     // debug
  //     std::cout << "    no potential next obs with high prob!!, using all" << std::endl;

  //     potential_next_obs = cells;
  //   }

  //   // debug
  //   // print out the next obs distribution
  //   std::cout << "  next obs distribution: " << std::endl;
  //   for( int i = 0; i < cells.size(); ++i ) {
  //     marked_grid_cell_t cell = cells[i];
  //     boost::optional<double> prob = next_observation_distribution( cell );
  //     if( prob ) {
  // 	std::cout << "    p( " << cell << " ) = " << *prob << std::endl;
  //     } else {
  // 	std::cout << "    p( " << cell << " ) = UNK! " << std::endl;
  //     }
  //   }
  //   std::cout.flush();

  //   // just pick the cell with the highest probability of being
  //   // in the shortest path
  //   double max_prob = *next_observation_distribution( potential_next_obs[0] );
  //   marked_grid_cell_t max_cell = potential_next_obs[0];
  //   for( size_t i = 0; i < potential_next_obs.size(); ++i ) {
  //     marked_grid_cell_t cell = potential_next_obs[i];
  //     double prob = *next_observation_distribution( cell );
  //     if( prob > max_prob ) {
  // 	max_prob = prob;
  // 	max_cell = cell;
  //     }
  //   }

  //   // debug
  //   std::cout << "    Picking next obs " << max_cell << ", p(..) = " << max_prob << std::endl;

  //   return max_cell;
  // }


  // //=======================================================================


  //=======================================================================

  marked_grid_cell_t 
  shortest_path_next_planner::
  choose_next_observation_cell()
  {

    scoped_context_switch context( chain_context( context_t( "greedy-shortest-path-choose-next-cell" ) ) );

    // estimate the most consistent workld in grids    
    std::vector<marked_grid_cell_t> binned_world =
      estimate_most_likely_world
      ( _point_process,
	_visited_grid,
	_sampler_planner_params.num_samples_of_point_sets,
	_sampler_planner_params.num_skip_between_point_set_samples);

    // Ok, now create virtual points at each cell
    std::vector<nd_point_t> virtual_points;
    for( size_t i = 0; i < binned_world.size(); ++i ) {
      nd_aabox_t region = _visited_grid.region( binned_world[i] );
      nd_point_t vp = region.start + 0.5 * ( region.end - region.start );
      virtual_points.push_back( vp );
    }
    
    // Now get next shortest path action in this virtual world
    boost::optional<marked_grid_cell_t> next_cell
      = next_observation_in_shortest_path( virtual_points,
					   _visited_grid,
					   _current_position );

    if( next_cell )
      return *next_cell;

    // Ok, hwat if the shortest path is invalid for some reason
    // Then jsut choose the first non-marked cell
    std::vector<marked_grid_cell_t> all_cells = _visited_grid.all_cells();
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      if( !_visited_grid( all_cells[i] ) ) {
	return all_cells[i];
      }
    }
    
    // getting here is unfortunate!
    std::cout << "AHHHHH! Why is there no unmarkes grid cell!" << std::endl;
    return _visited_grid.cell( _current_position );
  }

  //=======================================================================


  shortest_path_next_planner::
  shortest_path_next_planner
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const grid_planner_parameters_t& planner_params,
    const entropy_estimator_parameters_t& entropy_params,
    const sampler_planner_parameters_t& sampler_planner_params,
    const double& probability_next_action_is_in_shortest_path_threshold)
    : _point_process( process ),
      _planner_params( planner_params ),
      _entropy_params( entropy_params ),
      _sampler_planner_params( sampler_planner_params ),
      _probability_next_action_is_in_shortest_path_threshold(probability_next_action_is_in_shortest_path_threshold)
  {
    
    // create the grids with no marks using the process window
    _visited_grid = marked_grid_t<bool>( _point_process->window(),
					 _planner_params.grid_cell_size );
    _negative_observation_grid = marked_grid_t<bool>( _point_process->window(),
						      _planner_params.grid_cell_size );
    
  }

  //=======================================================================

  //=======================================================================
  
  std::string to_json( const grid_planner_parameters_t& p ) 
  {
    std::ostringstream oss;
    oss << "{ \"object_class\" : \"grid_planner_parameters_t\" , ";
    oss << "  \"burnin_mcmc_iterations\" : " << p.burnin_mcmc_iterations << " , ";
    oss << "  \"update_model_mcmc_iterations\" : " << p.update_model_mcmc_iterations << " , ";
    oss << "  \"grid_cell_size\" : " << p.grid_cell_size;
    oss << "}";
    return oss.str();
  }
  
  //=======================================================================

  std::string to_json( const sampler_planner_parameters_t& p )
  {
    std::ostringstream oss;
    oss << "{ \"object_class\" : \"sampler_planner_parameters_t\" , ";
    oss << "  \"num_samples_of_observations\" : " << p.num_samples_of_observations << " , ";
    oss << "  \"num_samples_of_point_sets\" : " << p.num_samples_of_point_sets;
    oss << "  \"num_skip_between_point_set_samples\" : " << p.num_skip_between_point_set_samples;
    oss << "}";
    return oss.str();
  }

  //=======================================================================

  std::string to_json( const entropy_estimator_parameters_t& p )
  {
    std::ostringstream oss;
    oss << "{ \"object_class\" : \"entropy_estimator_parameters_t\" , ";
    oss << "  \"num_samples\" : " << p.num_samples << " , ";
    oss << "  \"num_samples_to_skip\" : " << p.num_samples_to_skip << " , ";
    oss << "  \"histogram_grid_cell_size\" : " << p.histogram_grid_cell_size;
    oss << "}";
    return oss.str();
  }

  //=======================================================================

  std::string to_json( const marked_grid_cell_t& cell ) 
  {
    std::ostringstream oss;
    oss << "{ \"object_class\" : \"marked_grid_cell_t\" , ";
    oss << "  \"coordinate\" : [";
    for( size_t i = 0; i < cell.coordinate.size(); ++i ) {
      oss << cell.coordinate[i];
      if( i < cell.coordinate.size() - 1 ) {
	oss << ",";
      }
    }
    oss << "] }";
    return oss.str();
  }

  //=======================================================================

  std::string to_json( const marked_grid_t<bool>& grid )
  {
    std::ostringstream oss;
    oss << "{ \"object_class\" : \"marked_grid_t<bool>\" , ";
    oss << "  \"marked-cells\" : [ ";
    std::vector<marked_grid_cell_t> cells = grid.all_marked_cells();
    for( size_t i = 0; i < cells.size(); ++i ) {
      oss << to_json( cells[i] );
      if( i < cells.size() - 1 ) {
	oss << ",";
      }
    }
    oss << "] }";
    return oss.str();
  }

  //=======================================================================


  void 
  shortest_path_next_planner::print_shallow_trace( std::ostream& out ) const
  {

    out << "{ \"object_class\" : \"shortest_path_next_planner\" , ";
    out << "  \"current_position\" : " << to_json( _current_position ) << " , ";
    out << "  \"planner_params\" : " << to_json( _planner_params ) << " , ";
    out << "  \"entropy_params\" : " << to_json( _entropy_params ) << " , ";
    out << "  \"sample_planner_params\" : " << to_json( _sampler_planner_params ) << " , ";
    out << "  \"observations\" : [ ";
    for( size_t i = 0 ; i < _observations.size(); ++i ) {
      out << to_json( _observations[i]  );
      if( i < _observations.size() - 1 ) {
	out << ",";
      }
    }
    out << "], ";
    out << "  \"visited_grid\" : " << to_json( _visited_grid ) << " , ";
    out << "  \"negative_observation_grid\" : " << to_json( _negative_observation_grid ) << " , ";
    out << "  \"probability_next_action_is_in_shortest_path_threshold\" : " << _probability_next_action_is_in_shortest_path_threshold;
    out << "} ";
    
  }

  //=======================================================================
  //=======================================================================
  //=======================================================================

  // Description:
  // Crates a dataseries for a marked grid
  std::string
  create_dataseries( const marked_grid_t<double>& grid,
		     const ptree& extra_config,
		     const std::string& title,
		     const double& default_value = 0.0)
  {
    std::vector<data_point_t> data_points;
    for( auto cell : grid.all_cells_ordered() ) {
      double value = default_value;
      if( grid( cell ) ) {
	value = *grid( cell );
      }
      data_point_t dp = data_point_t( cell.coordinate[0],
				      cell.coordinate[1],
				      0.0 );
      dp.attributes.put( "value", value );
      data_points.push_back( dp );
    }

    ptree config = extra_config;
    std::string did
      = add_data_series( data_points,
			 config,
			 title );

    return did;
  }


  //=======================================================================
  
  std::string
  shortest_path_next_planner::plot( const std::string& title ) const
  {

    // get a plot for the model
    std::string model_plot = _point_process->plot( "model-for-" + title );
    
    // define a dataseries for the visited grid
    // std::string visited_ds 
    //   = create_dataseries( _visited_grid,
    // 			   ptree(),
    // 			   "visited-grid-of-" + title );
    
    // define a dataseries for the negative regions
    // std::string neg_ds
    //   = create_dataseries( _negative_observation_grid,
    // 			   ptree(),
    // 			   "neg-regions-of-" + title );
    
    // create a style for the visited and negative regions
    ptree visited_style;
    visited_style.put( "hidden", true );
    visited_style.put( "plot_prefix", "plot" );
    visited_style.put( "gnuplot.style", "image" );
    visited_style.add( "wanted_attributes.", "x" );
    visited_style.add( "wanted_attributes.", "y" );
    visited_style.add( "wanted_attributes.", "value" );
    ptree neg_style;
    neg_style.put( "hidden", true );
    neg_style.put( "plot_prefix", "plot" );
    neg_style.put( "gnuplot.style", "rgbalpha" );
    neg_style.add( "wanted_attributes.", "x" );
    neg_style.add( "wanted_attributes.", "y" );
    neg_style.add( "wanted_attributes.", "value" );
    neg_style.add( "wanted_attributes.", "zero" );
    neg_style.add( "wanted_attributes.", "zero" );
    neg_style.add( "wanted_attributes.", "value" );

    return model_plot;
  }


  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  


}
