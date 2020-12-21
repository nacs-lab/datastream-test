#

add_executable(cpu-single-thread-dry cpu-single-thread-dry.cpp)
target_link_libraries(cpu-single-thread-dry nacs-utils)

add_executable(cpu-single-thread-fill cpu-single-thread-fill.cpp)
target_link_libraries(cpu-single-thread-fill nacs-utils)
