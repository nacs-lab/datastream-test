#

add_definitions(-UNDEBUG)

add_executable(diff-uint64 diff-uint64.cpp)
target_link_libraries(diff-uint64 nacs-utils)

add_subdirectory(cpu)
if(ENABLE_OPENCL)
  add_subdirectory(opencl)
endif()
