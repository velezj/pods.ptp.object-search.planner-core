

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

using namespace math_core;
using namespace probability_core;
using namespace point_process_experiment_core;
using namespace ruler_point_process;
using namespace point_process_core;
using namespace planner_core;


void write_ruler_process_meta
( boost::shared_ptr<ruler_point_process_t>& process,
  std::ostream& fout_meta )
{
    // ok, writeo ut hte experiment details to the meta file
  fout_meta << "# OBSERVATIONS" << std::endl; 
 fout_meta << process->_state.observations.size() << std::endl;
  for( size_t i = 0; i < process->_state.observations.size(); ++i ) {
    fout_meta << process->_state.observations[i] << std::endl;
  }
  fout_meta << "# NEGATIVE OBSERVATION REGIONS" << std::endl;
  fout_meta << process->_state.negative_observations.size() << std::endl;
  for( size_t i = 0; i < process->_state.negative_observations.size(); ++i ) {
    fout_meta << process->_state.negative_observations[i] << std::endl;
  }
  fout_meta << "# INIT MIXTURE CLUSTERINGS" << std::endl;
  for( size_t i = 0; i < process->_state.observations.size(); ++i ) {
    fout_meta << process->_state.observation_to_mixture[i] << std::endl;
  }
  fout_meta << "# INIT MODEL" << std::endl;
  fout_meta << process->_state.model << std::endl;
  fout_meta << "# INIT MIXTURES" << std::endl;
  fout_meta << process->_state.mixture_gaussians.size() << std::endl;
  for( size_t i = 0; i < process->_state.mixture_gaussians.size(); ++i ) {
    fout_meta << "## MIXUTRE " << i << std::endl;
    fout_meta << "### SPREAD" << std::endl;
    fout_meta << process->_state.mixture_gaussians[i] << std::endl;
    fout_meta << "### PERIOD" << std::endl;
    fout_meta << process->_state.mixture_period_gammas[i] << std::endl;
    fout_meta << "### LENGTH" << std::endl;
    fout_meta << process->_state.mixture_ruler_length_gammas[i] << std::endl;
    fout_meta << "### START" << std::endl;
    fout_meta << process->_state.mixture_ruler_start_gaussians[i] << std::endl;
    fout_meta << "### DIRECTION" << std::endl;
    fout_meta << process->_state.mixture_ruler_direction_gaussians[i] << std::endl;
  }
  fout_meta << "# WINDOW" << std::endl;
  fout_meta << process->_state.window << std::endl;
  fout_meta << "# PRIORS" << std::endl;
  fout_meta << process->_state.model.prior_variance << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_p_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_p_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_q_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_q_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_r_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_r_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_s_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_period_s_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_p_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_p_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_q_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_q_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_r_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_r_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_s_shape << std::endl;
  fout_meta << process->_state.model.prior_ruler_length_s_rate << std::endl;
  fout_meta << process->_state.model.prior_ruler_start_mean << std::endl;
  fout_meta << process->_state.model.prior_ruler_start_variance << std::endl;
  fout_meta << process->_state.model.prior_ruler_direction_mean << std::endl;
  fout_meta << process->_state.model.prior_ruler_direction_variance << std::endl;
  fout_meta.flush();
}



std::vector<nd_aabox_t>
compute_empty_regions( const std::vector<nd_point_t>& points,
		       const nd_aabox_t& region,
		       const double epsilon = 1e-7)
{
  assert( region.start.n == 2 );
  assert( region.start.n > 0 );
  if( region.start.n != 2 ) {
    throw std::domain_error( "Canno compute empty regions, only implemented 2D case!" );
  }

  // ok, first we need to sort the points by their x coordinate
  // and another sorted by their y coordinates
  int dim = region.start.n;
  std::vector<nd_point_t> points_x( points );
  std::sort( points_x.begin(),
	     points_x.end(),
	     &point_compare_x );
  std::vector<nd_point_t> points_y( points );
  std::sort( points_y.begin(),
	     points_y.end(),
	     &point_compare_y );
  
  // Now, create the x and y lines
  std::vector< nd_point_t > x_line_ticks;
  std::vector< nd_point_t > y_line_ticks;
  x_line_ticks.push_back( region.start );
  y_line_ticks.push_back( region.start );
  for( size_t i = 0; i < points_x.size(); ++i ) {
    nd_vector_t x_margin = ( epsilon * vector( axis_direction( dim, 0 ) ) );
    x_line_ticks.push_back( points_x[i] + (-1.0 * x_margin) );
    x_line_ticks.push_back( points_x[i] + x_margin );
    nd_vector_t y_margin = ( epsilon * vector( axis_direction( dim, 1 ) ) );
    y_line_ticks.push_back( points_y[i] + (-1.0 * y_margin) );
    y_line_ticks.push_back( points_y[i] + y_margin );
  }
  x_line_ticks.push_back( region.end );
  y_line_ticks.push_back( region.end );
  
  // create all sub regions given the lines and return those which do not
  // have a point in them
  std::vector<nd_aabox_t> empty_regions;
  for( size_t xi = 0; xi < x_line_ticks.size() - 1; ++xi ) {
    nd_point_t xstart = x_line_ticks[xi];
    nd_point_t xend = x_line_ticks[xi+1];
    for( size_t yi = 0; yi < y_line_ticks.size() - 1; ++yi ) {
      nd_point_t ystart = y_line_ticks[yi];
      nd_point_t yend = y_line_ticks[yi+1];

      nd_point_t start = point( xstart.coordinate[0],
				ystart.coordinate[1] );
      nd_point_t end = point( xend.coordinate[0],
			      yend.coordinate[1] );
      nd_aabox_t box = aabox( start, end );
      bool empty = true;
      for( size_t k = 0; k < points.size(); ++k ) {
	if( is_inside( points[k], box ) ) {
	  empty = false;
	  break;
	}
      }
      if( empty ) {
	empty_regions.push_back( box );
      }
    }
  }
  
  return empty_regions;
}

// std::vector<nd_aabox_t> 
// compute_empty_regions( const std::vector<nd_point_t>& points,
// 		       const nd_aabox_t& region,
// 		       const double epsilon = 1e-7)
// {
//   assert( region.start.n == 1 );
//   if( region.start.n != 1 )
//     throw std::domain_error( "Cannot compute empty regions of non-1-dimensional regions" );


//   // add region start and end to poitns
//   std::vector<nd_point_t> p( points );
//   p.push_back( region.start );
//   p.push_back( region.end );

//   // first, sort points from low to high
//   std::sort( p.begin(), p.end(), &point_lexicographical_compare );
  
//   // ok, now create regions from start of region to first point,
//   // then first point to second point .. and so forth until region end
//   std::vector<nd_aabox_t> empty_regs;
//   for( size_t i = 0; i < p.size() - 1; ++i ) {
//     nd_aabox_t reg;
//     reg.n = region.n;

//     // the first point start is "flush" (so no epsilon offset)
//     if( i == 0 ) {
//       reg.start = p[i];
//     } else {
//       reg.start = p[i] + vector( epsilon );
//     }

//     // the last point end is "flush" (so no epsilon offset)
//     if( i == p.size() - 2 ) {
//       reg.end = p[i+1];
//     } else {
//       reg.end = p[i+1] + vector( -epsilon );
//     }

//     // add this empty space between points
//     empty_regs.push_back( reg );
//   }

//   // return emty space
//   return empty_regs;
// }



std::vector<nd_point_t>
observe_cell( shortest_path_next_planner& planner,
	      const marked_grid_cell_t& cell,
	      const std::vector<nd_point_t>& points )
{
  // Take any points inside the cell
  // and add as observations
  std::vector<nd_point_t> new_obs;
  nd_aabox_t region = planner.visisted_grid().region( cell );
  std::vector<nd_point_t> obs = planner.observations();
  for( std::size_t p_i = 0;
       p_i < points.size();
       ++p_i ) {
    
    // check if point is inside region and not already part of process
    nd_point_t point = points[ p_i ];
    if( is_inside( point, region ) &&
	std::find( obs.begin(),
		   obs.end(),
		   point ) == obs.end() ) {
      new_obs.push_back( point );
    }
  }
  
  // Ok, add new observation or a negative region if no new obs
  if( new_obs.empty() ) {
    planner.add_negative_observation( cell );
  } else {
    planner.add_observations( new_obs );
    
    // now add negative regions for the places in the cell without points
    std::vector<nd_aabox_t> empty_regs = compute_empty_regions( new_obs, region );
    for( size_t i = 0; i < empty_regs.size(); ++i ) {
      planner.add_empty_region( empty_regs[i] );
    }
  }
  
  return new_obs;
}


int main( int argc, char** argv )
{

  // get meta and trace files
  std::string dir = ".";
  if( argc > 1 ) {
    dir = argv[1];
  }
  std::ostringstream oss_meta;
  oss_meta << dir << "/" << "planner.meta";
  std::ofstream fout_meta( oss_meta.str().c_str() );
  std::ostringstream oss_trace;
  oss_trace << dir << "/" << "planner.trace";
  std::ofstream fout_trace( oss_trace.str().c_str() );
  
  // get seed if given
  unsigned int seed = 0;
  if( argc > 2 ) {
    std::istringstream iss_seed( argv[2] );
    iss_seed >> seed;
  }
  gsl_rng_set( global_rng(), seed );
  
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

  // some are known at start (cells 0 and 2 )
  std::vector<nd_point_t> init_points;
  init_points.push_back( points[0] );
  init_points.push_back( points[1] );
  init_points.push_back( points[5] );
  init_points.push_back( points[6] );
  
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
  				 init_points ) );
  boost::shared_ptr<mcmc_point_process_t> planner_process
    = boost::shared_ptr<mcmc_point_process_t>( process );
  

  // create a planner for it
  grid_planner_parameters_t planner_params;
  planner_params.burnin_mcmc_iterations = 100;
  planner_params.update_model_mcmc_iterations = 10;
  entropy_estimator_parameters_t entropy_params;
  entropy_params.num_samples = 10;
  sampler_planner_parameters_t sampler_planner_params;
  sampler_planner_params.num_samples_of_observations = 10;
  sampler_planner_params.num_samples_of_point_sets = 120;
  double prob_thresh = 0.6;
  shortest_path_next_planner 
    planner( planner_process,
	     planner_params,
	     entropy_params,
	     sampler_planner_params,
	     prob_thresh);

  // add observations initially to planner
  for( size_t i = 0 ; i < init_points.size(); ++i ) {
    planner._observations.push_back( init_points[i] );
  }

  // state for the planner controller
  std::size_t iteration = 0;
  std::vector<marked_grid_cell_t> chosen_cells;
  std::vector<nd_aabox_t> chosen_regions;
  std::vector<bool> chosen_region_negative;
  planner.set_current_position( window.start );

  // seed the initial seen cells
  chosen_cells.push_back( planner.visisted_grid().cell( point( 0.1, 0.1 ) ) );
  chosen_cells.push_back( planner.visisted_grid().cell( point( 2.1, 2.1 ) ) );
  chosen_regions.push_back( planner.visisted_grid().region( point( 0.1, 0.1 ) ) );
  chosen_regions.push_back( planner.visisted_grid().region( point( 2.1, 2.1 ) ) );
  chosen_region_negative.push_back( false );
  chosen_region_negative.push_back( false );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 0.1, 0.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 2.1, 2.1 ) ) );
  {
    std::vector<nd_point_t> p;
    p.push_back( init_points[0] );
    p.push_back( init_points[2] );
    nd_aabox_t region = planner.visisted_grid().region( point( 0.1, 0.1 ) );
    std::vector<nd_aabox_t> init_empty = compute_empty_regions( p, region );
    for( size_t i = 0; i < init_empty.size(); ++i ) {
      planner.add_empty_region( init_empty[i] );
    }
    p.clear();
    p.push_back( init_points[1] );
    p.push_back( init_points[3] );
    region = planner.visisted_grid().region( point( 2.1, 2.1 ) );
    init_empty = compute_empty_regions( p, region );
    for( size_t i = 0; i < init_empty.size(); ++i ) {
      planner.add_empty_region( init_empty[i] );
    }
  }
  planner.add_negative_observation( planner.visisted_grid().cell( point( 1.1, 1.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 3.1, 3.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 1.1, 1.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 3.1, 3.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 1.1, 3.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 3.1, 1.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 1.1, 3.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 3.1, 1.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 1.1, 2.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 3.1, 2.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 1.1, 2.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 3.1, 2.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 2.1, 1.1 ) ) );
  planner.add_negative_observation( planner.visisted_grid().cell( point( 2.1, 3.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 2.1, 1.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 2.1, 3.1 ) ) );
  
  
    

  // write the meta file
  fout_meta << "# RANDOM SEED" << std::endl;
  fout_meta << seed << std::endl;
  write_ruler_process_meta( process, fout_meta );
  fout_meta << "# PLANNER PARAMETERS" << std::endl;
  fout_meta << planner_params.burnin_mcmc_iterations << std::endl;
  fout_meta << planner_params.update_model_mcmc_iterations << std::endl;
  fout_meta << planner_params.grid_cell_size << std::endl;
  fout_meta << "# SAMPLER PLANNER PARAMETERS" << std::endl;
  fout_meta << sampler_planner_params.num_samples_of_observations << std::endl;
  fout_meta << sampler_planner_params.num_samples_of_point_sets << std::endl;
  fout_meta << "# ENTROPY PARAMETERS" << std::endl;
  fout_meta << entropy_params.num_samples << std::endl;
  fout_meta << entropy_params.num_samples_to_skip << std::endl;
  fout_meta << entropy_params.histogram_grid_cell_size << std::endl;
  fout_meta << "# INITIAL POINTS" << std::endl;
  for( size_t i = 0 ; i < init_points.size(); ++i ) {
    fout_meta << init_points[i] << std::endl;
  }
  fout_meta << "# INITIAL SEEN CELLS" << std::endl;
  fout_meta << planner.visisted_grid().cell( point( 0.1 ) ) << std::endl;
  fout_meta << planner.visisted_grid().cell( point( 2.1 ) ) << std::endl;

  // run the planner
  while( planner.observations().size() < points.size() ) {

    // Choose the next observation cell
    marked_grid_cell_t next_cell = 
      planner.choose_next_observation_cell();
    
    // Take any points inside the cell
    // and add as observations
    std::vector<nd_point_t> new_obs;
    nd_aabox_t region = planner.visisted_grid().region( next_cell );
    std::vector<nd_point_t> obs = planner.observations();
    for( std::size_t p_i = 0;
	 p_i < points.size();
	 ++p_i ) {

      // check if point is inside region and not already part of process
      nd_point_t point = points[ p_i ];
      if( is_inside( point, region ) &&
	  std::find( obs.begin(),
		     obs.end(),
		     point ) == obs.end() ) {
	new_obs.push_back( point );
      }
    }

    // update chosen cells and regions
    chosen_cells.push_back( next_cell );
    chosen_regions.push_back( region );
    
    // Ok, add new observation or a negative region if no new obs
    if( new_obs.empty() ) {
      planner.add_negative_observation( next_cell );
      chosen_region_negative.push_back( true );
    } else {
      
      // now add negative regions for the places in the cell without points
      std::vector<nd_aabox_t> empty_regs = compute_empty_regions( new_obs, region );
      for( size_t i = 0; i < empty_regs.size(); ++i ) {
	planner.add_empty_region( empty_regs[i] );
      }

      // and add teh actual observations 
      // (make sure this is AFTER the empty regions)
      planner.add_observations( new_obs );
      chosen_region_negative.push_back( false );      
    }

    // update position
    planner.set_current_position( region.start + (region.end - region.start) * 0.5 );


    // add the cell as visited to the planner
    planner.add_visited_cell( next_cell );

    // add to trace
    fout_trace << iteration << " "
	       << next_cell << " "
	       << new_obs.size() << " "
	       << planner.observations().size() << " "
	       << region << " ";
    for( size_t i = 0; i < new_obs.size(); ++i ) {
      fout_trace << new_obs[ i ] << " ";
    }
    fout_trace << std::endl;
    fout_trace.flush();

    // estiamte the entropy of the process after observations
    // double entropy =
    //   estimate_entropy_from_samples( entropy_params,
    // 				     process );

    // print the process model parameters to user
    std::cout << "[" << iteration << "]   " << std::endl;
    std::cout << process->_state << std::endl;
    
    // print status to user
    std::cout << "[" << iteration << "]   "
	      <<  "cell: " << next_cell 
	      << "  { #new= " << new_obs.size() << " total: " << planner.observations().size() << " }" 
      // << " entropy: " << entropy
	      << std::endl;
    std::cout << std::flush;
    
    // icrease iteration count
    ++iteration;
  }

  

  return 0;

}

