#

set(datastream_helpers_HDRS
  cpu_kernel.h
  test.h
  )
set(datastream_helpers_SRCS
  cpu_kernel.cpp
  test.cpp
  )
set(datastream_helpers_LINKS ${DEPS_LIBRARIES})

add_definitions("-\"DNACS_EXPORT_LIB_ds_helper()=\"")

add_library(datastream-helpers SHARED
  ${datastream_helpers_SRCS})
target_link_libraries(datastream-helpers PUBLIC ${datastream_helpers_LINKS})

set_target_properties(datastream-helpers PROPERTIES
  COMPILE_FLAGS "-fvisibility=hidden"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
