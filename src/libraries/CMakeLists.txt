set(MSD_LIBRARY_PRIMITIVES_SRC
  corefinementUtilities.cpp
  corefinementUtilities.h

  cylinderCorefinement.cpp
  cylinderCorefinement.h

  primitiveICP.cpp
  primitiveICP.h

  prismCorefinement.cpp
  prismCorefinement.h

  sphereCorefinement.cpp
  sphereCorefinement.h

  templatePrimitive.cpp
  templatePrimitive.h

  templateSpheres.cpp
  templateSpheres.h
)

set(MSD_LIBRARY_SKELETON_SRC
  skeletonExact.h
  skeletonExact.cpp

  medialSkeletalDiagramOpt.h
  medialSkeletalDiagramOpt.cpp
)

set(MSD_LIBRARY_SRC)

foreach(src ${MSD_LIBRARY_PRIMITIVES_SRC})
  list(APPEND MSD_LIBRARY_SRC "primitives/${src}")
endforeach()

foreach(src ${MSD_LIBRARY_SKELETON_SRC})
  list(APPEND MSD_LIBRARY_SRC "skeletons/${src}")
endforeach()

add_library(MSD_lib STATIC ${MSD_LIBRARY_SRC})
target_link_libraries(MSD_lib PUBLIC cgalInterface libiglInterface geogramInterface tetgenInterface configFileJSON nonlinearOptimization)
target_include_directories(MSD_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

source_group("primitives" FILES ${MSD_LIBRARY_PRIMITIVES_SRC})
source_group("skeletons" FILES ${MSD_LIBRARY_SKELETON_SRC})

set_property(TARGET MSD_lib PROPERTY FOLDER "projects/MSD")