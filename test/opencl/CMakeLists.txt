#

add_executable(opencl-dry-compute opencl-dry-compute.cpp)
target_link_libraries(opencl-dry-compute nacs-utils datastream-helpers)
