
#if !defined( __PLANNER_CORE_PLANNER_HPP__ )
#define __PLANNER_CORE_PLANNER_HPP__

#include <point-process-core/point_process.hpp>
#include <point-process-core/marked_grid.hpp>
#include <iosfwd>

namespace planner_core {

  using namespace point_process_core;


  // Description:
  // General parameters for a gird planner with a point process
  class grid_planner_parameters_t
  {
  public:

    std::size_t burnin_mcmc_iterations;
    std::size_t update_model_mcmc_iterations;

    double grid_cell_size;

    grid_planner_parameters_t()
      : burnin_mcmc_iterations( 100 ),
	update_model_mcmc_iterations( 10 ),
	grid_cell_size( 1.0 )
    {}

  };


  // Description:
  // General parameter for sampler planners (planners which must take future
  // observation samples)
  class sampler_planner_parameters_t
  {
  public:
    std::size_t num_samples_of_observations;
    std::size_t num_samples_of_point_sets;

    sampler_planner_parameters_t()
      : num_samples_of_observations( 100 ),
	num_samples_of_point_sets( 100 )
    {}
  };


  // Description:
  // A planner which works on a discrete marked grid of possible
  // observations to take 
  class grid_planner_t
  {
  public:

    // Description:
    // Returns the observation grid/ visited grid which has marked as true
    // all cells which have been observed
    virtual
    marked_grid_t<bool>
    visited_grid() const = 0;

    // Description:
    // Adds a given cell as already visisted
    virtual
    void add_visited_cell( const marked_grid_cell_t& cell ) = 0;
    
    // Description:
    // Add a negative region
    virtual
    void add_negative_observation( const marked_grid_cell_t& cell ) = 0;

    // Description:
    // Adds an "Empty" region, which is not a negative observation, rather
    // a region within an observed cell with points that has no points
    // hence the empty space of a cell with points
    virtual
    void add_empty_region( const nd_aabox_t& region ) = 0;

    // Description:
    // Add a set of observation poitns
    virtual
    void add_observations( const std::vector<math_core::nd_point_t>& obs ) = 0;

    // Description:
    // Returns the next observation to take
    virtual
    marked_grid_cell_t choose_next_observation_cell() = 0;

    // Description:
    // Sets the current position
    virtual
    void set_current_position( const math_core::nd_point_t& pos ) = 0;

    // Description:
    // Retruns the observations
    virtual
    std::vector<math_core::nd_point_t> observations() const = 0;

    // Description:
    // Returns true iff all cells have been visited
    virtual
    bool all_cells_visited() const = 0;

    // Description:
    // Prints out the "shallow" trace for this planner.
    // This is just a single-line (no newlione!_) string representation
    // of the parameters of this planner NOT including any deep objects
    // (hence the *shallow* part).  Only prints out things owned/operated by
    // this planner itself (so not the model, etc..)
    virtual
    void print_shallow_trace( std::ostream& out ) const = 0;

    // Descrioption:
    // Prints out a "shallow" trace for the model used by this planner.
    // This should just route to the model object's print_shallow_trace
    // method.
    virtual
    void print_model_shallow_trace( std::ostream& out ) const = 0;

  protected:

  };

}

#endif

