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
find_package(Doxygen)

# Options
option(FLT_BUILD_EXAMPLES "Build examples" ON)

# Core library
add_subdirectory(src)

# Examples
if(FLT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Documentation
if(DOXYGEN_FOUND)
    message("-- Doxygen found. Documentation will be built")
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    message("-- Configuring Doxyfile")
    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("-- Doxygen build started")

    add_custom_target(doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating Documentation with Doxygen"
        VERBATIM)
    install(DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/html DESTINATION share/doc)
else(DOXYGEN_FOUND)
    message("-- Doxygen not found. Documentation will not be generated.")
endif(DOXYGEN_FOUND)
