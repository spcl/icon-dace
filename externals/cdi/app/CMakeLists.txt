add_executable(cdi
  cdi.c printinfo.c
)

target_include_directories(cdi PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include
)

target_link_libraries(cdi PUBLIC cdilib)

target_compile_definitions(cdi PUBLIC HAVE_CONFIG_H)
