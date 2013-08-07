cmake_minimum_required(VERSION 2.6.0)

# pull in the pods macros. See cmake/pods.cmake for documentation
set(POD_NAME planner-core)
include(cmake/pods.cmake)

# automatically build LCM types.  This also defines a number of CMake
# variables, see cmake/lcmtypes.cmake for details
include(cmake/lcmtypes.cmake)
lcmtypes_build()


#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O2 -g")
#set(CMAKE_CXX_FLAGS "-std=c++0x -pedantic -Wall -O0 -g3")

add_definitions( -std=c++0x -pedantic -Wall )
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")


## the library
add_library( planner-core SHARED
  src/greedy_entropy_reduction_planner.cpp
  src/greedy_shortest_path_distribution_threshold_planner.cpp
  src/shortest_path_next_planner.cpp
)
pods_install_headers( 
  src/planner.hpp
  src/greedy_entropy_reduction_planner.hpp
  src/edmonds_optimum_branching.hpp
  src/edmonds_optimum_branching_impl.hpp
  src/greedy_shortest_path_distribution_threshold_planner.hpp
  src/shortest_path_next_planner.hpp
  DESTINATION
  planner-core
)
pods_use_pkg_config_packages(planner-core point-process-core lcmtypes_point-process-core boost-1.51.0)
pods_install_libraries( planner-core )
pods_install_pkg_config_file(planner-core
    CFLAGS
    LIBS -lplanner-core
    REQUIRES point-process-core boost-1.51.0
    VERSION 0.0.1)



add_subdirectory( test )
