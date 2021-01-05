cmake_minimum_required(VERSION 3.0)

project(datastream-test)

## Setup
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(GNUInstallDirs)

option(ENABLE_OPENCL "Enable OpenCL testing" On)

# Remove rdynamic
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

find_package(PkgConfig REQUIRED)
if(ENABLE_OPENCL)
  find_package(CLHPP)
  if(NOT CLHPP_FOUND)
    set(ENABLE_OPENCL Off FORCE)
  endif()
endif()
find_package(yaml-cpp REQUIRED)
pkg_check_modules(DEPS REQUIRED nacs-utils>=10.2)

add_definitions(-Wall -Wextra -g -pthread -fvisibility=hidden -Wconversion
  -fno-math-errno -Wno-sign-conversion -Wno-psabi -D_FILE_OFFSET_BITS=64
  ${DEPS_CFLAGS} ${TCDEP_CFLAGS})
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # We have GCC specific warning options that makes Clang unhappy.
  # Adding this unconditionally cause GCC to issue a warning if there are other warnings
  # so only add this conditionally for clang.
  add_definitions(-Wno-unknown-warning-option)
endif()
include_directories(${DEPS_INCLUDE_DIRS} ${YAML_CPP_INCLUDE_DIR})

set(LINKER_FLAGS "-Wl,--as-needed -Wl,--no-undefined -Wl,--gc-sections -pthread -fno-math-errno")
set(LINKER_FLAGS "${LINKER_FLAGS} -Wl,-Bsymbolic-function ${DEPS_LDFLAGS_OTHER}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17 -ggdb3 -ffp-contract=fast")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}")
set(CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LINKER_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(helpers)
add_subdirectory(test)
