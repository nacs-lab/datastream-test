#

add_executable(opencl-dry-compute opencl-dry-compute.cpp)
target_link_libraries(opencl-dry-compute nacs-utils datastream-helpers)

add_executable(opencl-compute-accuracy opencl-compute-accuracy.cpp)
target_link_libraries(opencl-compute-accuracy nacs-utils datastream-helpers -fopenmp)
target_compile_options(opencl-compute-accuracy PRIVATE -fopenmp -ffast-math)

add_executable(opencl-write-bandwidth opencl-write-bandwidth.cpp)
target_link_libraries(opencl-write-bandwidth nacs-utils datastream-helpers)

add_executable(opencl-write-bandwidth-mem-type opencl-write-bandwidth-mem-type.cpp)
target_link_libraries(opencl-write-bandwidth-mem-type nacs-utils datastream-helpers)

add_executable(opencl-parallel-rw opencl-parallel-rw.cpp)
target_link_libraries(opencl-parallel-rw nacs-utils datastream-helpers)

add_executable(opencl-rw-latency opencl-rw-latency.cpp)
target_link_libraries(opencl-rw-latency nacs-utils datastream-helpers)

add_executable(opencl-event-overhead opencl-event-overhead.cpp)
target_link_libraries(opencl-event-overhead nacs-utils datastream-helpers)
