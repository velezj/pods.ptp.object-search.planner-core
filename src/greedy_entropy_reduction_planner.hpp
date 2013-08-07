
#if !defined( __PLANNER_CORE_GREEDY_ENTROPY_REDUCTION_PLANNER_HPP__ )
#define __PLANNER_CORE_GREEDY_ENTROPY_REDUCTION_PLANNER_HPP__


#include "planner.hpp"
#include <point-process-core/entropy.hpp>

namespace planner_core {


  using namespace point_process_core;
  using namespace math_core;


  
  // Description:
  // A grid planner which tries to minimize the
  // entropy of hte point process model with the next observation
  // greedily (only looks at one observation in the future)
  class greedy_entropy_reduction_planner_t : public grid_planner_t
  {
  public:
    
    // Description:
    // Creates a new planner with the given point process as model
    // and entropy parameters
    greedy_entropy_reduction_planner_t
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const grid_planner_parameters_t& planner_params,
      const entropy_estimator_parameters_t& entropy_params,
      const sampler_planner_parameters_t& sampler_planner_params );
  

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
    // Sets the current position (currently IGNORES this)
    virtual
    void set_current_position( const nd_point_t& pos )
    {
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

  protected:

    
    // Description:
    // Estimate the entropy of the process after taking a future
    // observation region ( so the expected entropy after sampling
    // all points in hte given region as a grid cell)
    double estimate_entropy_after_future_observation_region
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const marked_grid_t<bool>& visited_grid,
      const marked_grid_cell_t& cell,
      const grid_planner_parameters_t& planner_params,
      const entropy_estimator_parameters_t& entropy_params,
      const sampler_planner_parameters_t& samples_planner_params ) const;

      
    
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

    
  };

}

#endif
