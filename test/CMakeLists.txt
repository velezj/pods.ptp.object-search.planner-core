
add_executable( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern.cpp)
pods_use_pkg_config_packages( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-igmm-point-process  lcmtypes_p2l-igmm-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern )



add_executable( planner-core-test-edmonds-branching
  test-edmonds-branching.cpp)
pods_use_pkg_config_packages( planner-core-test-edmonds-branching
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-igmm-point-process  lcmtypes_p2l-igmm-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  boost-1.54.0
  )
pods_install_executables( planner-core-test-edmonds-branching )





add_executable( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  test-greedy-shortest-path-distribution-threshold-planner.cpp)
pods_use_pkg_config_packages( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-igmm-point-process  lcmtypes_p2l-igmm-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( planner-core-test-greedy-shortest-path-distribution-threshold-planner )




add_executable( p2l-planner-core-experiment-1
  experiment-1.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-1
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-1 )

add_executable( p2l-planner-core-experiment-1c
  experiment-1c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-1c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-1c )


add_executable( p2l-planner-core-experiment-1-igmm
  experiment-1-igmm.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-1-igmm
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-igmm-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-1-igmm )




add_executable( p2l-planner-core-experiment-2
  experiment-2.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-2
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-2 )


add_executable( p2l-planner-core-experiment-2b
  experiment-2b.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-2b
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-2b )


add_executable( p2l-planner-core-experiment-2c
  experiment-2c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-2c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-2c )


add_executable( p2l-planner-core-experiment-2-igmm
  experiment-2-igmm.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-2-igmm
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-igmm-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-2-igmm )






add_executable( p2l-planner-core-experiment-3c
  experiment-3c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-3c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-3c )






add_executable( p2l-planner-core-experiment-4c
  experiment-4c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-4c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-4c )



add_executable( p2l-planner-core-experiment-5c
  experiment-5c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-5c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-5c )



add_executable( p2l-planner-core-test-console
  test-console.cpp)
pods_use_pkg_config_packages( p2l-planner-core-test-console
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-test-console )


add_executable( p2l-planner-core-experiment-6c
  experiment-6c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-6c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-6c )



add_executable( p2l-planner-core-experiment-7c
  experiment-7c.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-7c
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-7c )


add_executable( p2l-planner-core-experiment-8
  experiment-8.cpp)
pods_use_pkg_config_packages( p2l-planner-core-experiment-8
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core 
  p2l-probability-core lcmtypes_p2l-probability-core
  p2l-probability-core-graphics 
  p2l-ruler-point-process
  p2l-point-process-experiment-core
  p2l-planner-core
  )
pods_install_executables( p2l-planner-core-experiment-8 )
