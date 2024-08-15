set(MSD_LIBRARY_PRIMITIVES_SRC
  primitives/corefinementUtilities.cpp
  primitives/corefinementUtilities.h

  primitives/cylinderCorefinement.cpp
  primitives/cylinderCorefinement.h

  primitives/primitiveICP.cpp
  primitives/primitiveICP.h

  primitives/prismCorefinement.cpp
  primitives/prismCorefinement.h

  primitives/sphereCorefinement.cpp
  primitives/sphereCorefinement.h

  primitives/templatePrimitive.cpp
  primitives/templatePrimitive.h

  primitives/templateSpheres.cpp
  primitives/templateSpheres.h
)

set(MSD_LIBRARY_SKELETON_SRC
  skeletons/skeletonExact.h
  skeletons/skeletonExact.cpp

  skeletons/medialSkeletalDiagramOpt.h
  skeletons/medialSkeletalDiagramOpt.cpp

  skeletons/skeletonRVD.h
  skeletons/skeletonRVD.cpp
)

set(MSD_LIBRARY_SRC)

foreach(src ${MSD_LIBRARY_PRIMITIVES_SRC})
  list(APPEND MSD_LIBRARY_SRC "${src}")
endforeach()

foreach(src ${MSD_LIBRARY_SKELETON_SRC})
  list(APPEND MSD_LIBRARY_SRC "${src}")
endforeach()

add_library(MSD_lib STATIC ${MSD_LIBRARY_SRC})
target_link_libraries(MSD_lib PUBLIC cgalInterface libiglInterface geogramInterface tetgenInterface configFileJSON nonlinearOptimization nlopt)
target_include_directories(MSD_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

source_group("primitives" FILES ${MSD_LIBRARY_PRIMITIVES_SRC})
source_group("skeletons" FILES ${MSD_LIBRARY_SKELETON_SRC})

set_property(TARGET MSD_lib PROPERTY FOLDER "projects/MSD")