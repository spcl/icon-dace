add_custom_target("download_minimal_example_grids"
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/grid.1_0.nc ${CMAKE_CURRENT_BINARY_DIR}/grid.2_0.nc ${CMAKE_CURRENT_BINARY_DIR}/grid.2_1.nc
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_SOURCE_DIR}/download_minimal_example_grids.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
