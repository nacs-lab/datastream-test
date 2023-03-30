#

add_executable(cpu-single-thread-dry cpu-single-thread-dry.cpp)
target_link_libraries(cpu-single-thread-dry nacs-utils datastream-helpers)

add_executable(cpu-single-thread-fill cpu-single-thread-fill.cpp)
target_link_libraries(cpu-single-thread-fill nacs-utils datastream-helpers)

add_executable(cpu-write-bandwidth cpu-write-bandwidth.cpp)
target_link_libraries(cpu-write-bandwidth nacs-utils datastream-helpers)

add_executable(cpu-nt-write-bandwidth cpu-nt-write-bandwidth.cpp)
target_link_libraries(cpu-nt-write-bandwidth nacs-utils datastream-helpers)

add_executable(cpu-two-threads-ring-buffer cpu-two-threads-ring-buffer.cpp)
target_link_libraries(cpu-two-threads-ring-buffer nacs-utils datastream-helpers)

add_executable(cpu-two-threads-ring-buffer-compute cpu-two-threads-ring-buffer-compute.cpp)
target_link_libraries(cpu-two-threads-ring-buffer-compute nacs-utils datastream-helpers)

add_executable(cpu-multi-threads-compute cpu-multi-threads-compute.cpp)
target_link_libraries(cpu-multi-threads-compute nacs-utils datastream-helpers)

add_executable(cpu-compute-accuracy cpu-compute-accuracy.cpp)
target_link_libraries(cpu-compute-accuracy nacs-utils datastream-helpers -fopenmp)
target_compile_options(cpu-compute-accuracy PRIVATE -fopenmp -ffast-math)

add_executable(cpu-single-thread-full-sin cpu-single-thread-full-sin.cpp)
target_link_libraries(cpu-single-thread-full-sin nacs-utils datastream-helpers -fopenmp)
target_compile_options(cpu-single-thread-full-sin PRIVATE -fopenmp)
