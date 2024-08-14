set(COMPUTE_MSD_SRC computeMSD.cpp)

add_executable(computeMSD ${COMPUTE_MSD_SRC})
target_link_libraries(computeMSD PRIVATE MSD_lib)

set_property(TARGET computeMSD PROPERTY FOLDER "projects/MSD")
