cmake_minimum_required(VERSION 3.10)

project(FLT VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(MSVC)
    add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
    add_compile_options("-Wall" "-Wextra" "$<$<CONFIG:RELEASE>:-O3>")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options("-stdlib=libc++")
    else()
        # nothing special for gcc at the moment
    endif()
endif()

set(CMAKE_VERBOSE_MAKEFILE ON)

# 3rd party dependencies
find_package(embree 3.3.0 REQUIRED)

# Options
option(FLT_BUILD_EXAMPLES "Build examples" ON)

# Core library
add_subdirectory(src)

# Examples
if(FLT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()
