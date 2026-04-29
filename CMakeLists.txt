cmake_minimum_required(VERSION 4.0)
project(CFDSolvers C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(CFDSolvers src/main.c)
target_include_directories(CFDSolvers PRIVATE libraries/gsl/include)

if (WIN32)
    set(GSL_LOCAL_LIB_DIR "${CMAKE_SOURCE_DIR}/libraries/gsl/libwin")
else ()
    set(GSL_LOCAL_LIB_DIR "${CMAKE_SOURCE_DIR}/libraries/gsl/liblnx")
endif ()

find_library(GSL_LOCAL_LIBRARY NAMES gsl PATHS "${GSL_LOCAL_LIB_DIR}" NO_DEFAULT_PATH)
find_library(GSL_LOCAL_CBLAS_LIBRARY NAMES gslcblas PATHS "${GSL_LOCAL_LIB_DIR}" NO_DEFAULT_PATH)

if (GSL_LOCAL_LIBRARY AND GSL_LOCAL_CBLAS_LIBRARY)
    target_link_libraries(CFDSolvers PRIVATE "${GSL_LOCAL_LIBRARY}" "${GSL_LOCAL_CBLAS_LIBRARY}")
else ()
    find_package(GSL QUIET)
    if (GSL_FOUND)
        target_link_libraries(CFDSolvers PRIVATE GSL::gsl GSL::gslcblas)
    else ()
        message(WARNING "GSL headers are vendored in libraries/gsl/include, but no GSL binaries were found in ${GSL_LOCAL_LIB_DIR}. Build GSL into that folder or install GSL system-wide before linking.")
    endif ()
endif ()

if (NOT WIN32)
    target_link_libraries(CFDSolvers PRIVATE m)
endif ()
