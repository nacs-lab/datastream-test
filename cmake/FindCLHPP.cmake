#

find_package(OpenCL)

if(OpenCL_FOUND)
  find_path(CLHPP_INCLUDE_DIR
    NAMES
    CL/cl2.hpp
    PATHS ${OpenCL_INCLUDE_DIRS})
  set(CLHPP_INCLUDE_DIRS ${CLHPP_INCLUDE_DIR})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CLHPP
  FOUND_VAR CLHPP_FOUND
  REQUIRED_VARS OpenCL_FOUND CLHPP_INCLUDE_DIR)

mark_as_advanced(CLHPP_INCLUDE_DIR)

if(CLHPP_FOUND AND NOT TARGET OpenCL::CLHPP)
  add_library(OpenCL::CLHPP INTERFACE IMPORTED)
  target_link_libraries(OpenCL::CLHPP INTERFACE OpenCL::OpenCL)
  set_target_properties(OpenCL::CLHPP PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CLHPP_INCLUDE_DIRS}")
endif()
