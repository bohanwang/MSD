set(MSD_LIBRARY_PRIMITIVES_SRC
  corefinementUtilities.cpp
  corefinementUtilities.h

  cylinderCorefinement.cpp
  cylinderCorefinement.h

  sphereICP.cpp
  sphereICP.h

  primitiveICP.cpp
  primitiveICP.h

  prismCorefinement.cpp
  prismCorefinement.h

  sphereCorefinement.cpp
  sphereCorefinement.h

  sphereICPSolverBundle.cpp
  sphereICPSolverBundle.h

  templatePrimitive.cpp
  templatePrimitive.h

  
)

set(MSD_LIBRARY_SRC)

foreach(src ${MSD_LIBRARY_PRIMITIVES_SRC})
  list(APPEND MSD_LIBRARY_SRC "primitives/${src}")
endforeach()

add_library(MSD_lib STATIC ${MSD_LIBRARY_SRC})
target_link_libraries(MSD_lib PUBLIC cgalInterface libiglInterface nonlinearOptimization)


source_group("primitives" FILES ${MSD_LIBRARY_PRIMITIVES_SRC})