
#if !defined( __PLANNER_CORE_GREEDY_SHORTEST_PATH_DISTRIBUTION_THRESHOOLD_PLANNER_HPP__ )
#define __PLANNER_CORE_GREEDY_SHORTEST_PATH_DISTRIBUTION_THRESHOOLD_PLANNER_HPP__


#include "planner.hpp"
#include <point-process-core/entropy.hpp>


namespace planner_core {


  using namespace point_process_core;
  using namespace math_core;


  
  // Description:
  // A grid planner which tries to minimize the
  // shortest path given a shortest path distribution
  // with a particular threshold on the probability of the 
  // path being shroter than the taken action
  class greedy_shortest_path_distribution_threshold_planner 
    : public grid_planner_t
  {
  public:
    
    // Description:
    // Creates a new planner with the given point process as model
    // and entropy parameters
    greedy_shortest_path_distribution_threshold_planner
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const grid_planner_parameters_t& planner_params,
      const entropy_estimator_parameters_t& entropy_params,
      const sampler_planner_parameters_t& sampler_planner_params,
      const double& probability_next_action_is_in_shortest_path_threshold );
  

    // Description:
    // Returns the observation grid/ visited grid which has marked as true
    // all cells which have been observed
    virtual
    marked_grid_t<bool>
    visisted_grid() const
    {
      return _visited_grid;
    }

    // Description:
    // Adds a given cell as already visisted.
    virtual
    void add_visited_cell( const marked_grid_cell_t& cell )
    {
      _visited_grid.set( cell, true );
    }
    
    // Description:
    // Add a negative region
    virtual
    void add_negative_observation( const marked_grid_cell_t& cell )
    {
      _negative_observation_grid.set( cell, true );
      _point_process->add_negative_observation( _negative_observation_grid.region(cell) );
      _point_process->mcmc( _planner_params.update_model_mcmc_iterations );
    }

    // Description:
    // Adds an "Empty" region, which is not a negative observation, rather
    // a region within an observed cell with points that has no points
    // hence the empty space of a cell with points
    virtual
    void add_empty_region( const nd_aabox_t& region )
    {
      _point_process->add_negative_observation( region );
    }


    // Description:
    // Add a set of observation poitns
    virtual
    void add_observations( const std::vector<math_core::nd_point_t>& obs )
    {
      _point_process->add_observations( obs );
      _point_process->mcmc( _planner_params.update_model_mcmc_iterations );
      _observations.insert( _observations.end(), obs.begin(), obs.end() );
    }

    // Description:
    // Returns the next observation to take
    virtual
    marked_grid_cell_t choose_next_observation_cell();

    // Description:
    // Sets the current position
    virtual
    void set_current_position( const nd_point_t& pos )
    { 
      _current_position = pos;
    }

    // Description:
    // Retruns the observations
    virtual
    std::vector<math_core::nd_point_t> observations() const
    {
      return _observations;
    }

    // Description:
    // Returns true iff all cells have been visited
    virtual
    bool all_cells_visited() const
    {
      // TODO:
      // Write this!
      return false;
    }

  public:


    // Description:
    // Estiamte the distrubtion (by counts) of each potential next observation
    // cell given a particular new observatio nregion
    void
    estimate_shortest_path_next_observation_distribution_given_observation
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const marked_grid_t<bool>& visited_grid,
      const marked_grid_cell_t& cell,
      const grid_planner_parameters_t& planner_params,
      const sampler_planner_parameters_t& sampler_planner_params,
      const nd_point_t& current_position,
      marked_grid_t<double>& next_observation_distribution ) const;
    

    // Description:
    // Estimates the distribution of each potential observation cell
    // being the cell in the shortest path to observing all points.
    void
    estimate_shortest_path_next_observation_distribution
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const marked_grid_t<bool>& visisted_grid,
      const unsigned long num_samples_of_point_sets,
      const nd_point_t& current_position,
      marked_grid_t<double>& next_observation_distribution) const;


    // Description:
    // Returns the first observation in the 
    // shortest path sequence of observations to find
    // all of the points
    boost::optional<marked_grid_cell_t>
    next_observation_in_shortest_path
    ( const std::vector<nd_point_t>& points,
      const marked_grid_t<bool> visisted_grid,
      const nd_point_t& current_position) const;


    // Description:
    // the current position
    nd_point_t _current_position;
    
    // Description:
    // The point process used as model
    boost::shared_ptr<mcmc_point_process_t> _point_process;

    // Description:
    // The parameters for planning and entropy computations
    grid_planner_parameters_t _planner_params;
    entropy_estimator_parameters_t _entropy_params;
    sampler_planner_parameters_t _sampler_planner_params;

    // Description:
    // The observations so far
    std::vector<math_core::nd_point_t> _observations;
    
    // Description:
    // The grid of visited cells
    marked_grid_t<bool> _visited_grid;
        
    // Description:
    // The grid of negative observation
    marked_grid_t<bool> _negative_observation_grid;


    // Description:
    // Threshold of "good" potentail next actions in shortest path
    // as the probability that it is ion the shortest path
    double _probability_next_action_is_in_shortest_path_threshold;
    
  };

}

#endif
