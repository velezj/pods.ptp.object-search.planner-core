
#include "random_planner.hpp"
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include <math-core/utils.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <boost/pointer_cast.hpp>
#include <stdexcept>


#define VERBOSE false
#define PRINT_WARNINGS true

using namespace point_process_core;
using namespace math_core;

namespace planner_core {

  //=======================================================================
  //=======================================================================
  //=======================================================================


  //=======================================================================

  marked_grid_cell_t 
  random_planner_t::
  choose_next_observation_cell()
  {

    std::vector<marked_grid_cell_t> possible_cells;

    // Ok, what if the max cell  is invalid for some reason
    // Then jsut choose the first non-marked cell
    std::vector<marked_grid_cell_t> all_cells = _visited_grid.all_cells();
    for( size_t i = 0; i < all_cells.size(); ++i ) {
      if( !_visited_grid( all_cells[i] ) ) {
	possible_cells.push_back( all_cells[i] );
      }
    }

    // randomly choose a cell
    if( possible_cells.size() > 0 ) {
      std::random_shuffle( possible_cells.begin(), possible_cells.end() );
      return possible_cells[0];
    }
    
    // getting here is unfortunate!
    if( PRINT_WARNINGS ) {
      std::cout << "AHHHHH! Why is there no unmarkes grid cell!" << std::endl;
    }
    return _visited_grid.cell( _current_position );

  }

  //=======================================================================


  random_planner_t::
  random_planner_t
  ( const boost::shared_ptr<mcmc_point_process_t>& process,
    const grid_planner_parameters_t& planner_params )
    : _point_process( process ),
      _planner_params( planner_params )
  {
    
    // create the grids with no marks using the process window
    _visited_grid = marked_grid_t<bool>( _point_process->window(),
					 _planner_params.grid_cell_size );
    _negative_observation_grid = marked_grid_t<bool>( _point_process->window(),
						      _planner_params.grid_cell_size );
    
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
  


}
