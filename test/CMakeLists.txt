
add_executable( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern.cpp)
pods_use_pkg_config_packages( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  igmm-point-process  lcmtypes_igmm-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern )



add_executable( planner-core-test-edmonds-branching
  test-edmonds-branching.cpp)
pods_use_pkg_config_packages( planner-core-test-edmonds-branching
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  igmm-point-process  lcmtypes_igmm-point-process
  point-process-experiment-core
  planner-core
  boost-1.51.0
  )
pods_install_executables( planner-core-test-edmonds-branching )





add_executable( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  test-greedy-shortest-path-distribution-threshold-planner.cpp)
pods_use_pkg_config_packages( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  igmm-point-process  lcmtypes_igmm-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-test-greedy-shortest-path-distribution-threshold-planner )




add_executable( planner-core-experiment-1
  experiment-1.cpp)
pods_use_pkg_config_packages( planner-core-experiment-1
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-1 )

add_executable( planner-core-experiment-1c
  experiment-1c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-1c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-1c )


add_executable( planner-core-experiment-1-igmm
  experiment-1-igmm.cpp)
pods_use_pkg_config_packages( planner-core-experiment-1-igmm
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  igmm-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-1-igmm )




add_executable( planner-core-experiment-2
  experiment-2.cpp)
pods_use_pkg_config_packages( planner-core-experiment-2
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-2 )


add_executable( planner-core-experiment-2b
  experiment-2b.cpp)
pods_use_pkg_config_packages( planner-core-experiment-2b
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-2b )


add_executable( planner-core-experiment-2c
  experiment-2c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-2c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-2c )


add_executable( planner-core-experiment-2-igmm
  experiment-2-igmm.cpp)
pods_use_pkg_config_packages( planner-core-experiment-2-igmm
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  igmm-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-2-igmm )






add_executable( planner-core-experiment-3c
  experiment-3c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-3c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-3c )






add_executable( planner-core-experiment-4c
  experiment-4c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-4c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-4c )



add_executable( planner-core-experiment-5c
  experiment-5c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-5c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-5c )



add_executable( planner-core-test-console
  test-console.cpp)
pods_use_pkg_config_packages( planner-core-test-console
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-test-console )


add_executable( planner-core-experiment-6c
  experiment-6c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-6c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-6c )



add_executable( planner-core-experiment-7c
  experiment-7c.cpp)
pods_use_pkg_config_packages( planner-core-experiment-7c
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-7c )


add_executable( planner-core-experiment-8
  experiment-8.cpp)
pods_use_pkg_config_packages( planner-core-experiment-8
  math-core lcmtypes_math-core
  point-process-core lcmtypes_point-process-core
  probability-core lcmtypes_probability-core
  probability-core-graphics 
  ruler-point-process
  point-process-experiment-core
  planner-core
  )
pods_install_executables( planner-core-experiment-8 )
