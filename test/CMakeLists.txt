
add_executable( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern.cpp)
pods_use_pkg_config_packages( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process  lcmtypes_object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( test-greedy-entropy-reduction-planner-with-igmm-point-process-line-pattern )



add_executable( planner-core-test-edmonds-branching
  test-edmonds-branching.cpp)
pods_use_pkg_config_packages( planner-core-test-edmonds-branching
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process  lcmtypes_object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  boost-1.54.0
  )
pods_install_executables( planner-core-test-edmonds-branching )





add_executable( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  test-greedy-shortest-path-distribution-threshold-planner.cpp)
pods_use_pkg_config_packages( planner-core-test-greedy-shortest-path-distribution-threshold-planner
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process  lcmtypes_object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( planner-core-test-greedy-shortest-path-distribution-threshold-planner )




add_executable( object-search.planner-core-experiment-1
  experiment-1.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-1
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-1 )

add_executable( object-search.planner-core-experiment-1c
  experiment-1c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-1c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-1c )


add_executable( object-search.planner-core-experiment-1-igmm
  experiment-1-igmm.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-1-igmm
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-1-igmm )




add_executable( object-search.planner-core-experiment-2
  experiment-2.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-2
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-2 )


add_executable( object-search.planner-core-experiment-2b
  experiment-2b.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-2b
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-2b )


add_executable( object-search.planner-core-experiment-2c
  experiment-2c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-2c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-2c )


add_executable( object-search.planner-core-experiment-2-igmm
  experiment-2-igmm.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-2-igmm
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-2-igmm )






add_executable( object-search.planner-core-experiment-3c
  experiment-3c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-3c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-3c )






add_executable( object-search.planner-core-experiment-4c
  experiment-4c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-4c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-4c )



add_executable( object-search.planner-core-experiment-5c
  experiment-5c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-5c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-5c )



add_executable( object-search.planner-core-test-console
  test-console.cpp)
pods_use_pkg_config_packages( object-search.planner-core-test-console
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-test-console )


add_executable( object-search.planner-core-experiment-6c
  experiment-6c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-6c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-6c )



add_executable( object-search.planner-core-experiment-7c
  experiment-7c.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-7c
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-7c )


add_executable( object-search.planner-core-experiment-8
  experiment-8.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-8
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.ruler-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-8 )


add_executable( object-search.planner-core-experiment-8-igmm
  experiment-8-igmm.cpp)
pods_use_pkg_config_packages( object-search.planner-core-experiment-8-igmm
  object-search.math-core lcmtypes_object-search.math-core
  object-search.point-process-core 
  object-search.probability-core lcmtypes_object-search.probability-core
  object-search.probability-core-graphics 
  object-search.igmm-point-process
  object-search.point-process-experiment-core
  object-search.planner-core
  )
pods_install_executables( object-search.planner-core-experiment-8-igmm )
