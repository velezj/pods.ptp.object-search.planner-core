

#include <point-process-experiment-core/simulated_data.hpp>
#include <math-core/io.hpp>
#include <math-core/matrix.hpp>
#include <igmm-point-process/igmm_point_process.hpp>
#include <planner-core/greedy_shortest_path_distribution_threshold_planner.hpp>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <probability-core/distribution_utils.hpp>
#include <probability-core/core.hpp>
#include <limits>
#include <fstream>
#include <sstream>

using namespace math_core;
using namespace probability_core;
using namespace point_process_experiment_core;
using namespace igmm_point_process;
using namespace point_process_core;
using namespace planner_core;


void write_igmm_process_meta
( boost::shared_ptr<igmm_point_process_t>& process,
  const nd_aabox_t& window,
  std::ostream& fout_meta )
{
    // ok, writeo ut hte experiment details to the meta file
  fout_meta << "# OBSERVATIONS" << std::endl;
  fout_meta << process->_state.observations.size() << std::endl;
  for( int i = 0; i < process->_state.observations.size(); ++i ) {
    fout_meta << process->_state.observations[i] << std::endl;
  }
  fout_meta << "# NEGATIVE OBSERVATION REGIONS" << std::endl;
  fout_meta << process->_state.negative_observations.size() << std::endl;
  for( int i = 0; i < process->_state.negative_observations.size(); ++i ) {
    fout_meta << process->_state.negative_observations[i] << std::endl;
  }
  fout_meta << "# INIT MIXTURE CLUSTERINGS" << std::endl;
  for( int i = 0; i < process->_state.observations.size(); ++i ) {
    fout_meta << process->_state.observation_to_mixture[i] << std::endl;
  }
  fout_meta << "# INIT MODEL" << std::endl;
  fout_meta << process->_state.model << std::endl;
  fout_meta << "# INIT MIXTURES" << std::endl;
  fout_meta << process->_state.mixture_gaussians.size() << std::endl;
  for( int i = 0; i < process->_state.mixture_gaussians.size(); ++i ) {
    fout_meta << "## MIXUTRE " << i << std::endl;
    fout_meta << "### SPREAD" << std::endl;
    fout_meta << process->_state.mixture_gaussians[i] << std::endl;
    fout_meta << "### NUMBER" << std::endl;
    fout_meta << process->_state.mixture_poissons[i] << std::endl;
  }
  fout_meta << "# WINDOW" << std::endl;
  fout_meta << window << std::endl;
  fout_meta << "# PRIORS" << std::endl;
  fout_meta << process->_state.model.prior_mean << std::endl;
  fout_meta << process->_state.model.prior_variance << std::endl;
  fout_meta.flush();
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
  std::ofstream fout_meta( oss_meta.str() );
  std::ostringstream oss_trace;
  oss_trace << dir << "/" << "planner.trace";
  std::ofstream fout_trace( oss_trace.str() );

  // get seed if given
  unsigned int seed = 0;
  if( argc > 2 ) {
    std::istringstream iss_seed( argv[2] );
    iss_seed >> seed;
  }
  gsl_rng_set( global_rng(), seed );
  

  // craete the window (the range)
  nd_aabox_t window;
  window.n = 1;
  window.start = point( 0 );
  window.end = point( 10 );

  std::vector<nd_point_t> points;
  points.push_back( point( 0.1 ) );
  points.push_back( point( 3.1 ) );
  points.push_back( point( 6.1 ) );
  points.push_back( point( 9.1 ) );
  points.push_back( point( 0.9  ) );
  points.push_back( point( 3.9  ) );
  points.push_back( point( 6.9  ) );

  // some are known at start (cells 0 and 2 )
  std::vector<nd_point_t> init_points;
  init_points.push_back( points[0] );
  init_points.push_back( points[1] );
  init_points.push_back( points[4] );
  init_points.push_back( points[5] );
  
  // create the process
  int dim = 1;
  igmm_point_process_model_t model;
  model.alpha = 1;
  model.mean_distribution.dimension = 1;
  model.mean_distribution.means.push_back( 5 );
  model.mean_distribution.covariance = to_dense_mat( Eigen::MatrixXd::Identity(1,1) * 25.0 );
  model.precision_distribution.shape = 2;
  model.precision_distribution.rate = 0.25;
  model.num_points_per_gaussian_distribution.shape = 2;
  model.num_points_per_gaussian_distribution.rate = 0.25;
  model.prior_mean = 0;
  model.prior_variance = 1;
  boost::shared_ptr<igmm_point_process_t> process
    ( new igmm_point_process_t( window, model, init_points ));
  boost::shared_ptr<mcmc_point_process_t> planner_process
    = boost::shared_ptr<mcmc_point_process_t>( process );
  

  // create a planner for it
  grid_planner_parameters_t planner_params;
  planner_params.burnin_mcmc_iterations = 100;
  planner_params.update_model_mcmc_iterations = 30;
  entropy_estimator_parameters_t entropy_params;
  entropy_params.num_samples = 10;
  sampler_planner_parameters_t sampler_planner_params;
  sampler_planner_params.num_samples_of_observations = 10;
  sampler_planner_params.num_samples_of_point_sets = 10;
  double prob_thresh = 0.6;
  greedy_shortest_path_distribution_threshold_planner 
    planner( planner_process,
	     planner_params,
	     entropy_params,
	     sampler_planner_params,
	     prob_thresh);

  // add observations initially to planner
  for( int i = 0 ; i < init_points.size(); ++i ) {
    planner._observations.push_back( init_points[i] );
  }

  // state for the planner controller
  std::size_t iteration = 0;
  std::vector<marked_grid_cell_t> chosen_cells;
  std::vector<nd_aabox_t> chosen_regions;
  std::vector<bool> chosen_region_negative;
  planner.set_current_position( window.start );

  // seed the initial seen cells
  chosen_cells.push_back( planner.visisted_grid().cell( point( 0.1 ) ) );
  chosen_cells.push_back( planner.visisted_grid().cell( point( 3.1 ) ) );
  chosen_regions.push_back( planner.visisted_grid().region( point( 0.1 ) ) );
  chosen_regions.push_back( planner.visisted_grid().region( point( 3.1 ) ) );
  chosen_region_negative.push_back( false );
  chosen_region_negative.push_back( false );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 0.1 ) ) );
  planner.add_visited_cell( planner.visisted_grid().cell( point( 3.1 ) ) );

  // write the meta file
  fout_meta << "# RANDOM SEED" << std::endl;
  fout_meta << seed << std::endl;
  write_igmm_process_meta( process, window, fout_meta );
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
  for( int i = 0 ; i < init_points.size(); ++i ) {
    fout_meta << init_points[i] << std::endl;
  }
  fout_meta << "# INITIAL SEEN CELLS" << std::endl;
  fout_meta << planner.visisted_grid().cell( point( 0.1 ) ) << std::endl;
  fout_meta << planner.visisted_grid().cell( point( 3.1 ) ) << std::endl;

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
    for( int i = 0; i < new_obs.size(); ++i ) {
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

