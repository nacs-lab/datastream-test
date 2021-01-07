#

add_executable(opencl-dry-compute opencl-dry-compute.cpp)
target_link_libraries(opencl-dry-compute nacs-utils datastream-helpers)

add_executable(opencl-compute-accuracy opencl-compute-accuracy.cpp)
target_link_libraries(opencl-compute-accuracy nacs-utils datastream-helpers -fopenmp)
target_compile_options(opencl-compute-accuracy PRIVATE -fopenmp -ffast-math)
