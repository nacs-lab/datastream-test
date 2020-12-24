#

add_executable(cpu-single-thread-dry cpu-single-thread-dry.cpp)
target_link_libraries(cpu-single-thread-dry nacs-utils datastream-helpers)

add_executable(cpu-single-thread-fill cpu-single-thread-fill.cpp)
target_link_libraries(cpu-single-thread-fill nacs-utils datastream-helpers)

add_executable(cpu-write-bandwidth cpu-write-bandwidth.cpp)
target_link_libraries(cpu-write-bandwidth nacs-utils datastream-helpers)

add_executable(cpu-nt-write-bandwidth cpu-nt-write-bandwidth.cpp)
target_link_libraries(cpu-nt-write-bandwidth nacs-utils datastream-helpers)
