

#include <point-process-experiment-core/simulated_data.hpp>
#include <math-core/io.hpp>
#include <math-core/matrix.hpp>
#include <igmm-point-process/igmm_point_process.hpp>
#include <planner-core/greedy_shortest_path_distribution_threshold_planner.hpp>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <probability-core/distribution_utils.hpp>
#include <probability-core-graphics/lcmgl_distributions.hpp>
#include <limits>

using namespace math_core;
using namespace probability_core;
using namespace probability_core::graphics;
using namespace point_process_experiment_core;
using namespace igmm_point_process;
using namespace point_process_core;
using namespace planner_core;


void colorize( long i, double& r, double& g, double& b )
{
  static long color_n = 11;
  static double color_r[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  static double color_g[] = { 0.3, 0.2, 0.1, 0.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4 };
  static double color_b[] = { 0.2, 0.3, 0.5, 0.7, 0.2, 0.3, 0.5, 0.7, 0.2, 0.3, 0.5 };
  
  // colorize negatiuve values 
  if( i < 0 ) {
    colorize( std::numeric_limits<long>::max() + i, r, g, b );
    return;
  }
  
  // for positive indices
  r = color_r[ i % color_n ];
  g = color_g[ i % color_n ];
  b = color_b[ i % color_n ];
}


void draw_igmm( bot_lcmgl_t* lcmgl,
		boost::shared_ptr<igmm_point_process_t>& proccess )
{
  igmm_point_process_state_t state = proccess->_state;
  
  // draw the mixtures
  double r,g,b;
  double zero_p[] = { 0,0,0 };
  double x_scale = 1;
  double y_scale = 3;
  for( std::size_t mi = 0; mi < state.mixture_gaussians.size(); ++mi ) {
    colorize( mi, r,g,b );
    lcmglPushMatrix();
    lcmglScalef( x_scale, y_scale, 1.0 );
    lcmglColor3f( r, g, b );
    draw_distribution( lcmgl, state.mixture_gaussians[mi] );
    lcmglColor3f( 1, 1, 1 );
    bot_lcmgl_text( lcmgl, zero_p, "Gauss" );
    lcmglPopMatrix();
    lcmglPushMatrix();
    lcmglTranslated( 0.0, y_scale, 0.0 );
    lcmglScalef( x_scale, y_scale, 1.0 );
    lcmglTranslated( state.mixture_gaussians[mi].means[0], 0.0, 0.0 );
    lcmglColor3f( r, g, b );
    draw_distribution( lcmgl, state.mixture_poissons[mi] );
    lcmglColor3f( 1,1,1 );
    bot_lcmgl_text( lcmgl, zero_p, "Poss" );
    lcmglPopMatrix();
  }
  lcmglPushMatrix();
  lcmglTranslated( 0, 2 * y_scale, 0 );
  lcmglScalef( x_scale, y_scale, 1 );
  lcmglColor3f( 0,0,0 );
  draw_distribution( lcmgl, state.model.mean_distribution );
  lcmglColor3f( 1,1,1 );
  bot_lcmgl_text( lcmgl, zero_p, "Mean" );
  lcmglPopMatrix();
  lcmglTranslated( 0, 3 * y_scale, 0 );
  lcmglScalef( x_scale, y_scale, 1 );
  lcmglColor3f( 0,0,0 );
  draw_distribution( lcmgl, state.model.precision_distribution );
  lcmglColor3f( 1,1,1 );
  bot_lcmgl_text( lcmgl, zero_p, "Prec" );
  lcmglPopMatrix();
  lcmglTranslated( 0, 4 * y_scale, 0 );
  lcmglScalef( x_scale, y_scale, 1 );
  lcmglColor3f( 0,0,0 );
  draw_distribution( lcmgl, state.model.num_points_per_gaussian_distribution );
  lcmglColor3f( 1,1,1 );
  bot_lcmgl_text( lcmgl, zero_p, "Num" );
  lcmglPopMatrix();
  
}


void draw_cells( bot_lcmgl_t* lcmgl,
		 const greedy_shortest_path_distribution_threshold_planner& planner,
		 const std::vector<marked_grid_cell_t>& cells,
		 const std::vector<nd_aabox_t>& regions,
		 const std::vector<bool>& region_neg )
{
  double y_extra = 0.2;
  char buf[512];
  for( size_t i = 0; i < cells.size(); ++i ) {
    sprintf( buf, "%zi", i );
    nd_aabox_t reg = regions[i];
    bool neg = region_neg[i];
    if( neg ) {
      lcmglColor3f( 1, 0, 0.3 );
    } else {
      lcmglColor3f( 0, 1, 0.3 );
    }
    lcmglBegin( LCMGL_QUADS );
    lcmglVertex3d( reg.start.coordinate[0], -y_extra, 0.1 );
    lcmglVertex3d( reg.end.coordinate[0]  , -y_extra, 0.1 );
    lcmglVertex3d( reg.end.coordinate[0]  ,  y_extra, 0.1 );
    lcmglVertex3d( reg.start.coordinate[0],  y_extra, 0.1 );
    lcmglEnd();
    double x_mid = reg.start.coordinate[0] + (reg.end.coordinate[0] - reg.start.coordinate[0]) / 2;
    double loc[] = { x_mid, -y_extra - 1, 0.1 };
    bot_lcmgl_text( lcmgl, loc, buf );
  }
}


void draw_points( bot_lcmgl_t* lcmgl,
		  const std::vector<nd_point_t>& points,
		  const double& radius = 0.1,
		  const double& height = 0.5)
{
  for( size_t i = 0; i < points.size(); ++i ) {
    double loc[] = { points[i].coordinate[0], 0, height }; 
    lcmglDisk( loc, 0, radius );
  }
}

int main( int argc, char** argv )
{

  lcm_t* lcm = lcm_create(NULL);
  bot_lcmgl_t* igmm_lcmgl = bot_lcmgl_init( lcm, "IGMM" );
  bot_lcmgl_t* planner_lcmgl = bot_lcmgl_init( lcm, "PLANNER" );
  bot_lcmgl_t* points_lcmgl = bot_lcmgl_init( lcm, "TRUE-POINTS" );
  bot_lcmgl_t* lcmgl;


  // craete the window (the range)
  nd_aabox_t window;
  window.n = 1;
  window.start = point( 0 );
  window.end = point( 10 );

  // simulate some points
  std::vector<nd_point_t> points 
    = simulate_line_point_clusters_gaussian_spread_poisson_size
    ( window,
      3,
      1.5,
      3.0 );

  // draw the tru points
  lcmgl = points_lcmgl;
  lcmglColor3f( 0,0,0  );
  draw_points( lcmgl, points );
  bot_lcmgl_switch_buffer( points_lcmgl );

  
  // now create a new point process;
  igmm_point_process_model_t model;
  model.alpha = 1;
  model.mean_distribution.dimension = 1;
  model.mean_distribution.means.push_back( 5 );
  model.mean_distribution.covariance = to_dense_mat( Eigen::MatrixXd::Identity(1,1) * 1.0 );
  model.precision_distribution.shape = 2;
  model.precision_distribution.rate = 0.25;
  model.num_points_per_gaussian_distribution.shape = 2;
  model.num_points_per_gaussian_distribution.rate = 0.25;
  model.prior_mean = 0;
  model.prior_variance = 1;
  boost::shared_ptr<igmm_point_process_t> igmm_process( new igmm_point_process_t( window, model, std::vector<nd_point_t>() ));
  boost::shared_ptr<mcmc_point_process_t> process( igmm_process );
  

  // create a planner for it
  grid_planner_parameters_t planner_params;
  entropy_estimator_parameters_t entropy_params;
  sampler_planner_parameters_t sampler_planner_params;
  double prob_thresh = 0.6;
  greedy_shortest_path_distribution_threshold_planner 
    planner( process,
	     planner_params,
	     entropy_params,
	     sampler_planner_params,
	     prob_thresh);

  // tell the user how many poitns we generates
  std::cout << "Total Points: " << points.size() << std::endl;

  // run the planner
  std::size_t iteration = 0;
  std::vector<marked_grid_cell_t> chosen_cells;
  std::vector<nd_aabox_t> chosen_regions;
  std::vector<bool> chosen_region_negative;
  planner.set_current_position( window.start );
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

    // draw the igmm model
    draw_igmm( igmm_lcmgl, igmm_process );

    // draw the cells
    draw_cells( planner_lcmgl,
		planner,
		chosen_cells,
		chosen_regions,
		chosen_region_negative );

    // swap all lcmgl buffers
    bot_lcmgl_switch_buffer( igmm_lcmgl );
    bot_lcmgl_switch_buffer( planner_lcmgl );

    // add the cell as visited to the planner
    planner.add_visited_cell( next_cell );

    // estiamte the entropy of the process after observations
    double entropy =
      estimate_entropy_from_samples( entropy_params,
				     process );

    // print the process model parameters to user
    std::cout << "[" << iteration << "]   "
	      << "alpha: " << igmm_process->_state.model.alpha 
	      << " mu: " << igmm_process->_state.model.mean_distribution
	      << " prec: " << igmm_process->_state.model.precision_distribution
	      << " num: " << igmm_process->_state.model.num_points_per_gaussian_distribution << std::endl;
    
    // print status to user
    std::cout << "[" << iteration << "]   "
	      <<  "cell: " << next_cell 
	      << "  { #new= " << new_obs.size() << " total: " << planner.observations().size() << " }" 
	      << " entropy: " << entropy
	      << std::endl;
    std::cout << std::flush;
    
    // icrease iteration count
    ++iteration;
  }

  

  return 0;

}

