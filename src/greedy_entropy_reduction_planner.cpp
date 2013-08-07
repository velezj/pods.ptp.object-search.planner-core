
#include "greedy_entropy_reduction_planner.hpp"
#include <math-core/utils.hpp>
#include <math-core/geom.hpp>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math-core/io.hpp>


using namespace math_core;
using namespace point_process_core;


namespace planner_core {


  //=========================================================================

  greedy_entropy_reduction_planner_t::greedy_entropy_reduction_planner_t
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const grid_planner_parameters_t& planner_params,
    const entropy_estimator_parameters_t& entropy_params,
    const sampler_planner_parameters_t& sampler_planner_params)
    : _point_process( process ),
      _planner_params( planner_params ),
      _entropy_params( entropy_params ),
      _sampler_planner_params( sampler_planner_params )
  {
    
    // create the grids with no marks using the process window
    _visited_grid = marked_grid_t<bool>( _point_process->window(),
					 _planner_params.grid_cell_size );
    _negative_observation_grid = marked_grid_t<bool>( _point_process->window(),
						      _planner_params.grid_cell_size );
    
  }

  //=========================================================================

  marked_grid_cell_t
  greedy_entropy_reduction_planner_t::choose_next_observation_cell()
  {

    // debug
    //std::cout << "  choosing next cell..." << std::endl << std::flush;
    
    // For every cell which is not visited, compute the entropy
    // of taking an observation there
    // and keep track of lowest
    std::vector<std::pair<marked_grid_cell_t,double> > entropies;
    std::vector<marked_grid_cell_t> all_cells 
      = _visited_grid.all_cells();
    std::size_t best_index = 0;
    for( std::size_t cell_i = 0; cell_i < all_cells.size(); ++cell_i ) {
      marked_grid_cell_t cell = all_cells[ cell_i ];

      // skip visisted cells
      if( _visited_grid( cell ) )
	continue;
      
      // enstimate the entropy of observing the cell
      double entropy = 
	estimate_entropy_after_future_observation_region
	( _point_process,
	  _visited_grid,
	  cell,
	  _planner_params,
	  _entropy_params,
	  _sampler_planner_params);
      
      entropies.push_back( std::make_pair( cell, entropy ) );

      // debug
      //std::cout << "  estiamted entropy (" << cell << ") = " << entropy << std::endl << std::flush;

      // keep track of min entropy
      if( best_index != cell_i ) { // test for first case
	if( entropy < entropies[ best_index ].second ) {
	  best_index = entropies.size() - 1; // the current index
	}
      }
    }

    // debug
    std::cout << "  estimated entropy after obs " 
	      << entropies[ best_index ].first << " : "
	      << entropies[ best_index ].second << std::endl;

    // ok, just return the cell for the best index
    return entropies[ best_index ].first;
  }

  //=========================================================================

  double 
  greedy_entropy_reduction_planner_t::
  estimate_entropy_after_future_observation_region
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const marked_grid_t<bool>& visited_grid,
    const marked_grid_cell_t& cell,
    const grid_planner_parameters_t& planner_params,
    const entropy_estimator_parameters_t& entropy_params,
    const sampler_planner_parameters_t& sampler_planner_params ) const
  {

    // debug
    //std::cout << "    future obs: " << cell << std::endl << std::flush;
    
    // First, clone the process
    boost::shared_ptr<mcmc_point_process_t> process_clone = process->clone();

    // get the know observations of hte process
    std::vector<nd_point_t> original_observations = process_clone->observations();
    
    // Ok, get the region for hte cell
    nd_aabox_t region = visited_grid.region( cell );

    // We are going to sample some observations and store hte entoyp of 
    // hte process if those observations were taken
    std::vector<double> entropies;
    
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

      // debug
      //std::cout << "      [" << sample_obs_i << "] mixed with " << new_obs.size() << " new obs" << std::endl << std::flush;
      
      // Store the entropy of this new process state
      double entropy = 
	estimate_entropy_from_samples( entropy_params, 
 				       process_with_new_observations );
      entropies.push_back( entropy );

      // debug
      //std::cout << "    entropy: " << entropy << std::endl << std::flush;
	
    }


    // Return the expected entropy from samples
    return mean( entropies );
    
  }


  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================
  //=========================================================================


}
