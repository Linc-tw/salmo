

    ######################################################
    ##  CMakeLists.txt for Salmo v0.0                   ##
    ##  Version 2020.07.09                              ##
    ##                                                  ##
    ##  Copyright (C) 2020 - Chieh-An Lin               ##
    ##  GNU GPLv3 - https://www.gnu.org/licenses/       ##
    ######################################################


cmake_minimum_required(VERSION 2.6)
project(salmo)

option(CPP "Compile with Linc's path" OFF)

## Use pkg-config to find modules
find_package(PkgConfig REQUIRED)
pkg_check_modules(REQ_PKG REQUIRED gsl cfitsio chealpix)
#pkg_check_modules(MPI mpi)

## Link the required modules
include_directories(${REQ_PKG_INCLUDE_DIRS})
link_directories(${REQ_PKG_LIBRARY_DIRS})

## Optional module healpix_cxx, link if found
if(CPP)
  pkg_check_modules(HEALPIX_CXX healpix_cxx>=3.30)
  if(HEALPIX_CXX_FOUND)
    add_definitions(-D __MFP_USE_HEALPIX_CXX__)
    include_directories(${HEALPIX_CXX_INCLUDE_DIRS})
    link_directories(${HEALPIX_CXX_LIBRARY_DIRS})
  endif()
endif()

## Optional module mpi, link if found
#if(MPI_FOUND)
#  add_definitions(-D __MFP_USE_MPI__)
#  include_directories(${MPI_INCLUDE_DIRS})
#  link_directories(${MPI_LIBRARY_DIRS})
#endif()

## Define cflags
set(DEF_C_FLAGS "-Wall -Wpedantic -Wno-unused-result -Wno-maybe-uninitialized -O3 -ffast-math")
set(ADD_C_FLAGS "-Wno-unused-variable")
set(CMAKE_C_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")
if(CPP)
  set(CMAKE_CXX_FLAGS "${DEF_C_FLAGS} ${ADD_C_FLAGS}")
endif()

## Build executable
if(CPP)
  file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/source/*.c" "${PROJECT_SOURCE_DIR}/source/*.cpp")
  add_executable(salmo ${SRC_FILES})
  target_link_libraries(salmo
    ${REQ_PKG_LIBRARIES}
    ${HEALPIX_CXX_LIBRARIES}
  #  ${MPI_LIBRARIES}
  )
else()
  file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/source/*.c")
  add_executable(salmo ${SRC_FILES})
  target_link_libraries(salmo
    ${REQ_PKG_LIBRARIES}
    ${HEALPIX_CXX_LIBRARIES}
  )
endif()

