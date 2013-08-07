

#include <point-process-experiment-core/simulated_data.hpp>
#include <math-core/io.hpp>
#include <math-core/matrix.hpp>
#include <ruler-point-process/ruler_point_process.hpp>
#include <planner-core/shortest_path_next_planner.hpp>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <probability-core/distribution_utils.hpp>
#include <probability-core/core.hpp>
#include <limits>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <math-core/utils.hpp>



using namespace math_core;
using namespace probability_core;
using namespace point_process_experiment_core;
using namespace ruler_point_process;
using namespace point_process_core;
using namespace planner_core;


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
  for( int i = 0; i < all_cells.size(); ++i ) {

    // debug 
    //std::cout << "cell = " << all_cells[i] << std::endl;

    marked_grid_cell_t cell = all_cells[i];
    bool marked = (std::find( full_cells.begin(),
			      full_cells.end(),
			      cell ) != full_cells.end());
    
    // write out empoty of filled cell
    if( marked ) {
      std::cout << "v";
    } else {
      std::cout << "-";
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
      
      // debug
      //std::cout << "  [" << k << "] diff: " << diff << " grid_diff: " << grid_diff << " cs: " << cell_size << " (" << start_point << "," << cell_point << ")" << std::endl;
      
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
	
	// debug
	//std::cout << "  [" << k << "] dir: " << dir << ", nc: " << (cell_point + dir) << std::endl;
      } 
    }

    // debug
    //std::cout << "bc: " << barrier_count << " b: " << barrier_idx << " ec: " << end_count << " e: " << end_idx << std::endl;
    
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



int main()
{


  // craete the window (the range)
  nd_aabox_t window;
  window.n = 2;
  window.start = point( 0.0, 0.0 );
  window.end = point( 10.0, 10.0 );

  // create some points
  std::vector<nd_point_t> points;
  points.push_back( point( 0.25, 0.25 ) );
  points.push_back( point( 2.25, 2.25 ) );
  points.push_back( point( 4.25, 4.25 ) );
  points.push_back( point( 6.25, 6.25 ) );
  points.push_back( point( 8.25, 8.25 ) );
  points.push_back( point( 0.55, 0.55 ) );
  points.push_back( point( 2.55, 2.55 ) );
  points.push_back( point( 4.55, 4.55 ) );
  points.push_back( point( 6.55, 6.55 ) );
  points.push_back( point( 8.55, 8.55 ) );

    // create the process
  int dim = 2;
  ruler_point_process_model_t model;
  model.prior_ruler_start_mean = zero_point( dim );
  model.prior_ruler_direction_mean = zero_point( dim );
  model.alpha = 1;
  model.precision_distribution.shape = 500000;
  model.precision_distribution.rate = 1000;
  model.period_distribution.p = pow(2.15,40);
  model.period_distribution.q = 2.15*40;
  model.period_distribution.r = 40;
  model.period_distribution.s = 40;
  model.ruler_length_distribution.p = pow(10,40);
  model.ruler_length_distribution.q = 10 * 40;
  model.ruler_length_distribution.r = 40;
  model.ruler_length_distribution.s = 40;
  model.ruler_start_mean_distribution.dimension = dim;
  model.ruler_start_mean_distribution.means.push_back( 0.4 );
  model.ruler_start_mean_distribution.means.push_back( 0.4 );
  model.ruler_start_mean_distribution.covariance = to_dense_mat( Eigen::MatrixXd::Identity(dim,dim) * (0.4*0.4) );
  model.ruler_start_precision_distribution.shape = 5000;
  model.ruler_start_precision_distribution.rate = 100;
  model.ruler_direction_mean_distribution.dimension = dim;
  model.ruler_direction_mean_distribution.means.push_back( 10 );
  model.ruler_direction_mean_distribution.means.push_back( 10 );
  model.ruler_direction_mean_distribution.covariance = to_dense_mat( Eigen::MatrixXd::Identity(dim,dim) * 1 );
  model.ruler_direction_precision_distribution.shape = 5000;
  model.ruler_direction_precision_distribution.rate = 100;
  boost::shared_ptr<ruler_point_process_t> process = 
    boost::shared_ptr<ruler_point_process_t>
    ( new ruler_point_process_t( window,
  				 model,
  				 points ) );
  boost::shared_ptr<mcmc_point_process_t> planner_process
    = boost::shared_ptr<mcmc_point_process_t>( process );

  // create a planner for it
  grid_planner_parameters_t planner_params;
  planner_params.burnin_mcmc_iterations = 1;
  planner_params.update_model_mcmc_iterations = 1;
  entropy_estimator_parameters_t entropy_params;
  entropy_params.num_samples = 10;
  sampler_planner_parameters_t sampler_planner_params;
  sampler_planner_params.num_samples_of_observations = 10;
  sampler_planner_params.num_samples_of_point_sets = 5;
  double prob_thresh = 0.6;
  shortest_path_next_planner 
    planner( planner_process,
	     planner_params,
	     entropy_params,
	     sampler_planner_params,
	     prob_thresh);

  
  // get the visited grid
  marked_grid_t<bool> visited_grid = planner.visisted_grid();
  
  // print out to console
  spnp_console_print_point_set( points, visited_grid );

  return 0;
}
