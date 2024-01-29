# expose functions for downstream projects
file(COPY Scripts/ComInAddTest.cmake
  DESTINATION ${CMAKE_BINARY_DIR})
install(FILES Scripts/ComInAddTest.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")
