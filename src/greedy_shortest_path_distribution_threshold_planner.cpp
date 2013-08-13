
#include "greedy_shortest_path_distribution_threshold_planner.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include "edmonds_optimum_branching.hpp"
#include <point-process-core/context.hpp>
#include <iostream>
#include <iomanip>
#include <algorithm>


using namespace point_process_core;

namespace planner_core {


  //=======================================================================

  void console_print_point_set( const std::vector<nd_point_t>& points,
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

    // go from cells to single number (hack)
    std::vector<int> marks;
    for( size_t i = 0; i < full_cells.size(); ++i ) {
      marks.push_back( (int)(full_cells[i].coordinate[0]) );
    }

    // sort marks
    std::sort( marks.begin(), marks.end() );
    
    // "print" out the range (1D)
    for( int i = 0; i < 11; ++i ) {
      if( std::find( marks.begin(), marks.end(), i ) == marks.end() ) {
	std::cout << "_";
      } else {
	std::cout << "v";
      }
      if( i == 5 ) {
	std::cout << " | ";
      }
    }
    std::cout << std::endl;
    std::cout.flush();
  }

  //=======================================================================

  void console_print_next_obs_dist( const marked_grid_t<double>& grid )
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

  boost::optional<marked_grid_cell_t>
  greedy_shortest_path_distribution_threshold_planner
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
  greedy_shortest_path_distribution_threshold_planner::
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
	= process->sample();

      // skip updating distribution if no points in sample
      if( sample_point_set.empty() ) {

	// debug 
	std::cout << " #0 ";
	std::cout.flush();

	continue;
      }

      // print this sample to the console
      console_print_point_set( sample_point_set, visited_grid );
      
      // Calcualte the next observation in the shortest path thorugh
      // points
      boost::optional<marked_grid_cell_t> next_obs 
	= next_observation_in_shortest_path( sample_point_set,
					     visited_grid,
					     current_position );
      
      // skip updating distribution if no next observation found
      if( !next_obs ) {
	continue;
      }
      
      // Ok, increment the mark at this cell
      boost::optional<double> m = next_observation_distribution( *next_obs );
      if( m ) {
	next_observation_distribution.set( *next_obs, *m + 1.0 );
      } else {
	next_observation_distribution.set( *next_obs, 1.0 );
      }
    }

    // // Ok, now divide all counts by number of samples
    // std::vector<marked_grid_cell_t> cells = next_observation_distribution.all_marked_cells();
    // for( size_t i = 0; i < cells.size(); ++i ) {
    //   next_observation_distribution.set( cells[i],
    // 					 *next_observation_distribution( cells[i] ) / (double)num_samples_of_point_sets );
    // }

  }


  //=======================================================================

  marked_grid_cell_t 
  greedy_shortest_path_distribution_threshold_planner::
  choose_next_observation_cell()
  {

    scoped_context_switch context( chain_context( context_t( "greedy-shortest-path-choose-next-cell" ) ) );
    
    // For all non-visited cells,
    // compute the distribution over the expected
    // next action for shortest path to all points
    // given we take the observation
    marked_grid_t<double> next_observation_distribution = _visited_grid.copy_structure<double>();
    std::vector<marked_grid_cell_t> all_cells 
      = _visited_grid.all_cells();
    for( std::size_t cell_i = 0; cell_i < all_cells.size(); ++cell_i ) {
      marked_grid_cell_t cell = all_cells[ cell_i ];
      
      // skip visisted cells
      if( _visited_grid( cell ) )
	continue;
      
      
      // estimate the next observation distribution given an observation
      estimate_shortest_path_next_observation_distribution_given_observation
	( _point_process,
	  _visited_grid,
	  cell,
	  _planner_params,
	  _sampler_planner_params, 
	  _current_position,
	  next_observation_distribution );
    }

    // normalize the cells to be probabilities
    all_cells = next_observation_distribution.all_marked_cells();
    double sum = 0;
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      sum += *next_observation_distribution(all_cells[i]);
    }
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      next_observation_distribution.set( all_cells[i], *next_observation_distribution(all_cells[i]) / sum );
    }
    
    
    // Ok, now get the set of actions which have a high enough
    // probability
    std::vector<marked_grid_cell_t> cells = next_observation_distribution.all_marked_cells();
    std::vector<marked_grid_cell_t> potential_next_obs;
    potential_next_obs = cells;

    // debug
    // print out the next obs distribution
    std::cout << "  next obs distribution: " << std::endl;
    for( size_t i = 0; i < cells.size(); ++i ) {
      marked_grid_cell_t cell = cells[i];
      boost::optional<double> prob = next_observation_distribution( cell );
      if( prob ) {
	std::cout << "    p( " << cell << " ) = " << *prob << std::endl;
      } else {
	std::cout << "    p( " << cell << " ) = UNK! " << std::endl;
      }
    }
    std::cout.flush();

    // just pick the cell with the highest probability of being
    // in the shortest path
    double max_prob = *next_observation_distribution( potential_next_obs[0] );
    marked_grid_cell_t max_cell = potential_next_obs[0];
    size_t max_cell_i = 0;
    while( _visited_grid( max_cell ) && (max_cell_i+1) < potential_next_obs.size() ) {
      ++max_cell_i;
      max_prob = *next_observation_distribution( potential_next_obs[max_cell_i] );
      max_cell = potential_next_obs[max_cell_i];
    }
    for( size_t i = 0; i < potential_next_obs.size(); ++i ) {
      marked_grid_cell_t cell = potential_next_obs[i];
      double prob = *next_observation_distribution( cell );
      if( prob > max_prob && !_visited_grid( cell ) ) {
	max_prob = prob;
	max_cell = cell;
      }
    }

    // debug
    std::cout << "    Picking next obs " << max_cell << ", p(..) = " << max_prob << std::endl;

    return max_cell;
  }


  //=======================================================================

  greedy_shortest_path_distribution_threshold_planner::
  greedy_shortest_path_distribution_threshold_planner
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

  void
  greedy_shortest_path_distribution_threshold_planner::
  estimate_shortest_path_next_observation_distribution_given_observation
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const marked_grid_t<bool>& visited_grid,
    const marked_grid_cell_t& cell,
    const grid_planner_parameters_t& planner_params,
    const sampler_planner_parameters_t& sampler_planner_params,
    const nd_point_t& current_position,
    marked_grid_t<double>& next_observation_distribution ) const
  {

    scoped_context_switch context( chain_context( context_t( "next-obs-dist-given-obs" ) ) );
    
    // First, clone the process
    boost::shared_ptr<mcmc_point_process_t> process_clone = process->clone();
    
    // get the know observations of hte process
    std::vector<nd_point_t> original_observations = process_clone->observations();
    
    // Ok, get the region for hte cell
    nd_aabox_t region = visited_grid.region( cell );

    // make copy of visited grid and "visit" the observedcell
    marked_grid_t<bool> visited_grid_copy = visited_grid;
    visited_grid_copy.set( cell, true );

    // debug
    std::cout << "  estimating sp next obs, future obs cell: " << cell << region << std::endl;
    std::cout.flush();

    // Sample some points sets from the process
    // And make a copy of the process seeing those points which lie in 
    // the region, then compute the entropy of the resulting process state.
    for( std::size_t sample_obs_i = 0; sample_obs_i < sampler_planner_params.num_samples_of_observations; ++sample_obs_i ) {

      // sample from the process
      std::vector<nd_point_t> point_set = 
	process_clone->sample_and_step();

      // Now keep those points which are in the observation
      // and new
      std::vector<nd_point_t> new_obs;
      for( std::size_t p_i = 0; p_i < point_set.size(); ++p_i ) {
	nd_point_t point = point_set[ p_i ];

	// if this is a new point within region, set as new obs
	if( is_inside( point, region ) &&
	    std::find( original_observations.begin(),
		       original_observations.end(),
		       point ) == original_observations.end() ) {
	  new_obs.push_back( point );
	}
      }


      // if there are no new obs, add a negative observation otherwise
      // add the new observations to a clone
      boost::shared_ptr<mcmc_point_process_t> process_with_new_observations
	= process_clone->clone();
      if( new_obs.empty() ) {
	process_with_new_observations->add_negative_observation( region );
      } else {
	process_with_new_observations->add_observations( new_obs );
      }
      
      // MCMC step the process with the new observations so that it mixes
      process_with_new_observations->mcmc( planner_params.update_model_mcmc_iterations );

      // add coutns to the distribution of next shortest path observations
      estimate_shortest_path_next_observation_distribution
	( process_with_new_observations,
	  visited_grid_copy,
	  sampler_planner_params.num_samples_of_point_sets,
	  current_position,
	  next_observation_distribution );

      // debug
      //std::cout << ".";
      std::cout << "    ---- " << sample_obs_i <<  std::endl;
      console_print_next_obs_dist( next_observation_distribution );
      std::cout << "    ---- " << sample_obs_i <<  std::endl;
      std::cout.flush();
    }

    // debug
    std::cout << std::endl;
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
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  //=======================================================================
  


}
