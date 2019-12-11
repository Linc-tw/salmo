

    ######################################################
    ##  CMakeLists.txt for mockFootprint v0.0		##
    ##  Version 2019.01.31				##
    ##							##
    ##  Copyright (C) 2019 - Chieh-An Lin		##
    ##  GNU GPLv3 - https://www.gnu.org/licenses/	##
    ######################################################


cmake_minimum_required(VERSION 2.6)
project(mockFootprint)

## Use pkg-config to find modules
find_package(PkgConfig REQUIRED)
pkg_check_modules(REQ_PKG REQUIRED gsl fftw3 cfitsio chealpix)
pkg_check_modules(HEALPIX_CXX healpix_cxx)
#pkg_check_modules(MPI mpi)

## Link the required modules
include_directories(${REQ_PKG_INCLUDE_DIRS})
link_directories(${REQ_PKG_LIBRARY_DIRS})

## Optional module healpix_cxx, link if found
if (HEALPIX_CXX_FOUND)
  add_definitions(-D __MFP_USE_HEALPIX_CXX__)
  include_directories(${HEALPIX_CXX_INCLUDE_DIRS})
  link_directories(${HEALPIX_CXX_LIBRARY_DIRS})
endif()

## Optional module mpi, link if found
#if (MPI_FOUND)
#  add_definitions(-D __MFP_USE_MPI__)
#  include_directories(${MPI_INCLUDE_DIRS})
#  link_directories(${MPI_LIBRARY_DIRS})
#endif()

## Define cflags
set(DEF_C_FLAGS "-Wall -Wpedantic -Wno-unused-result -Wno-maybe-uninitialized -O3 -ffast-math")
set(ADD_C_FLAGS "-Wno-unused-variable")
set(CMAKE_C_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")

## Build executable
file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/source/*.c" "${PROJECT_SOURCE_DIR}/source/*.cpp")
add_executable(mockFootprint ${SRC_FILES})
target_link_libraries(mockFootprint
  ${REQ_PKG_LIBRARIES}
  ${HEALPIX_CXX_LIBRARIES}
#  ${MPI_LIBRARIES}
)

