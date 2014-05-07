cmake_minimum_required(VERSION 2.6.0)

# pull in the pods macros. See cmake/pods.cmake for documentation
set(POD_NAME p2l-planner-core)
include(cmake/pods.cmake)

# automatically build LCM types.  This also defines a number of CMake
# variables, see cmake/lcmtypes.cmake for details
include(cmake/lcmtypes.cmake)
lcmtypes_build()


#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O2 -g")
#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O0 -g3")

add_definitions( -std=c++0x -Wall -fdiagnostics-show-option -Wno-unused-local-typedefs -fPIC )
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")

option ( USE_PEDANTIC "Turn on -pendantic mode in gcc. This will spit out *lots* of warnings from lcm :-(, but hopefully none from the rest of the code" OFF)
if( USE_PEDANTIC )
    add_definitions( -pedantic )
endif (USE_PEDANTIC)


## the library
add_library( p2l-planner-core SHARED
  src/greedy_entropy_reduction_planner.cpp
  src/greedy_shortest_path_distribution_threshold_planner.cpp
  src/shortest_path_next_planner.cpp
  src/one_action_entropy_reduction_planner.cpp
)
pods_install_headers( 
  src/planner.hpp
  src/greedy_entropy_reduction_planner.hpp
  src/edmonds_optimum_branching.hpp
  src/edmonds_optimum_branching_impl.hpp
  src/greedy_shortest_path_distribution_threshold_planner.hpp
  src/shortest_path_next_planner.hpp
  src/coverage_planner.hpp
  src/one_action_entropy_reduction_planner.hpp
  DESTINATION
  planner-core
)
pods_use_pkg_config_packages(p2l-planner-core 
    p2l-point-process-core 
    boost-1.54.0
    plot-server-api-internal
    plot-server-api
    plot-server-util
    p2l-ruler-point-process
    p2l-igmm-point-process
)
pods_install_libraries( p2l-planner-core )
pods_install_pkg_config_file(p2l-planner-core
    CFLAGS
    LIBS -lp2l-planner-core
    REQUIRES p2l-point-process-core boost-1.54.0 plot-server-api-internal plot-server-api plot-server-util p2l-ruler-point-process p2l-igmm-point-process
    VERSION 0.0.1)



#add_subdirectory( test )
