# To use CUDA by enabling language
cmake_minimum_required(VERSION 3.18)
project(NAUNET VERSION 0.3.0 LANGUAGES CXX)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_definitions(-DNAUNET_DEBUG)
endif()

message(STATUS "CMAKE_SYSTEM_NAME          = ${CMAKE_SYSTEM_NAME}"         )
message(STATUS "CMAKE_CXX_COMPILER_ID      = ${CMAKE_CXX_COMPILER_ID}"     )
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_COMPILER_IS_GNUCXX   = ${CMAKE_COMPILER_IS_GNUCXX}"  )
message(STATUS "CMAKE_BUILD_TYPE           = ${CMAKE_BUILD_TYPE}"          )

set(SUNDIALS_MIN_VERSION "5.6.1")
set(SUNDIALS_DIR "" CACHE PATH "Path to SUNDIALSConfig.cmake")
find_package(SUNDIALS ${SUNDIALS_MIN_VERSION} REQUIRED CONFIG)

message(STATUS "SUNDIALS_FOUND             = ${SUNDIALS_FOUND}")

list(APPEND SUNDIALSLIB
    SUNDIALS::cvode
    SUNDIALS::sunlinsolklu
    SUNDIALS::sunmatrixsparse
    )

# find_package(CUDA REQUIRED)
# find_package(CUDAToolkit REQUIRED)

# if(CMAKE_BUILD_TYPE STREQUAL Debug)
#     cmake_print_properties(TARGETS SUNDIALS::cvode_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
#     cmake_print_properties(TARGETS SUNDIALS::nvecserial_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
# endif()

set(NVIDIA_GPU_COUNT 0)


message(STATUS "NVIDIA_GPU_COUNT           = " ${NVIDIA_GPU_COUNT})

set(MAKE_STATIC ON CACHE BOOL "")
set(MAKE_SHARED ON CACHE BOOL "")
set(MAKE_PYTHON OFF CACHE BOOL "")
set(MAKE_TEST ON CACHE BOOL "")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(MAKE_PYTHON)
    set(PYMODNAME "pynaunet" CACHE STRING "")
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
endif()

add_subdirectory(src)

if(MAKE_TEST)
    enable_testing()
    add_subdirectory(tests)
endif()
