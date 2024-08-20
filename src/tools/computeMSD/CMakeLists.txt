set(COMPUTE_MSD_SRC computeMSD.cpp)

add_executable(computeMSD ${COMPUTE_MSD_SRC})
target_link_libraries(computeMSD PRIVATE MSD_lib)
target_link_libraries(computeMSD PRIVATE argparse::argparse)

set_property(TARGET computeMSD PROPERTY FOLDER "projects/MSD")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set_target_properties(computeMSD PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/bin)
  set_target_properties(computeMSD PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/bin)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set_target_properties(computeMSD PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${PGO_BINARY_ROOT}/bin/Debug)
  set_target_properties(computeMSD PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PGO_BINARY_ROOT}/bin/Release)
endif()
