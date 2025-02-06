# expose functions for downstream projects
file(COPY Scripts/ComInAddTest.cmake
  DESTINATION ${PROJECT_BINARY_DIR})
install(FILES Scripts/ComInAddTest.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")
file(COPY Scripts/predef_replay_data.cmake
  DESTINATION ${PROJECT_BINARY_DIR})
install(FILES Scripts/predef_replay_data.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}/cmake")
