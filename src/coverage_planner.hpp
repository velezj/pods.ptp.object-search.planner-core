
#if !defined( __P2L_PLANNER_CORE_coverage_planner_HPP__ )
#define __P2L_PLANNER_CORE_coverage_planner_HPP__


#include "planner.hpp"

namespace planner_core {


  using namespace point_process_core;
  using namespace math_core;
  
  
  
  // Description:
  // a gridplanner which does simple coverage
  class coverage_planner 
    : public grid_planner_t
  {
  public:
    
    // Description:
    // Creates a new planner with the given grid
    coverage_planner
    ( const boost::shared_ptr<mcmc_point_process_t>& process,
      const grid_planner_parameters_t& planner_params )
      : _point_process( process ), 
	_planner_params( planner_params ),
	_current_position( zero_point(2) )
    {
      _visited_grid = marked_grid_t<bool>( _point_process->window(),
					   _planner_params.grid_cell_size );
      _negative_observation_grid = marked_grid_t<bool>( _point_process->window(),
							_planner_params.grid_cell_size );
    }

    virtual ~coverage_planner() { };
  

    // Description:
    // Returns the observation grid/ visited grid which has marked as true
    // all cells which have been observed
    virtual
    marked_grid_t<bool>
    visited_grid() const
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
    }

    // Description:
    // Adds an "Empty" region, which is not a negative observation, rather
    // a region within an observed cell with points that has no points
    // hence the empty space of a cell with points
    virtual
    void add_empty_region( const nd_aabox_t& region )
    {
    }


    // Description:
    // Add a set of observation poitns
    virtual
    void add_observations( const std::vector<math_core::nd_point_t>& obs )
    {
      _observations.insert( _observations.end(), obs.begin(), obs.end() );
    }

    // Description:
    // Returns the next observation to take
    virtual
    marked_grid_cell_t choose_next_observation_cell()
    {
      for( auto cell : _visited_grid.all_cells_ordered() ) {
	if( !_visited_grid( cell ) ) {
	  return cell;
	}
      }
      std::cout << "AHHH, there are NO MORE CELLS!" << std::endl;
      return _visited_grid.cell( _current_position );
    }

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

    // Description:
    // Print outs a shallow trace for this planner
    virtual
    void print_shallow_trace( std::ostream& out ) const
    {
    }

    // Description:
    // print shallow trace of the model
    virtual
    void print_model_shallow_trace( std::ostream& out ) const
    {
    }

    // Description:
    // Get/Set the grid_planner_parameters_t parameters
    virtual
    grid_planner_parameters_t get_grid_planner_parameters() const
    {
      return _planner_params;
    }
    virtual
    void set_grid_planner_parameters( const grid_planner_parameters_t& p )
    {
      _planner_params = p;
    }

  public:


    // Description:
    // the current position
    nd_point_t _current_position;
    
    // Description:
    // The point process used as model
    boost::shared_ptr<mcmc_point_process_t> _point_process;

    // Description:
    // The parameters for planning and entropy computations
    grid_planner_parameters_t _planner_params;

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
