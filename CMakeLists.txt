cmake_minimum_required(VERSION 3.12)

project(FLT VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Options
option(FLT_BUILD_EXAMPLES "Build examples" ON)
option(FLT_BUILD_TESTS "Build tests" ON)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executable")
set(INSTALL_DOC_DIR share/doc CACHE PATH "Installation directory for documentation")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/flt)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

if(MSVC)
    add_compile_options("/W4" "$<$<CONFIG:RELEASE>:/O2>")
else()
    add_compile_options("-Wall" "-Wextra" "$<$<CONFIG:RELEASE>:-O3>")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options("-stdlib=libc++")
        add_link_options("-stdlib=libc++")
    else()
        # nothing special for gcc at the moment
    endif()
endif()

# set(CMAKE_VERBOSE_MAKEFILE ON)

# 3rd party dependencies
find_package(embree REQUIRED)
# Set Embree MAJOR version
set(FLT_COMPILE_DEFINITIONS "EMBREE_VERSION=${EMBREE_VERSION_MAJOR}")
add_compile_definitions(${FLT_COMPILE_DEFINITIONS})

find_package(OpenMP REQUIRED)
# Optional for documentation generation
find_package(Doxygen)

# Core library
add_subdirectory(src)

# Examples
if(FLT_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Tests
if(FLT_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
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

    add_custom_command(OUTPUT doc_doxygen/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating Documentation with Doxygen"
        VERBATIM)
    add_custom_target(doc ALL DEPENDS doc_doxygen/html/index.html)
    install(DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/doc_doxygen/html DESTINATION share/doc)
else(DOXYGEN_FOUND)
    message("-- Doxygen not found. Documentation will not be generated.")
endif(DOXYGEN_FOUND)

# -- Export configuration
# Straight from How to create a ProjectConfig.cmake file

# Add all targets to the build-tree export set
export(TARGETS flt flt_static
  FILE "${PROJECT_BINARY_DIR}/FLTTargets.cmake")

# Export the package for use from the build-tree.
export(PACKAGE FLT)

# Create the FLTConfig.cmake and FLTConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(FLTConfig.cmake.in
  "${PROJECT_BINARY_DIR}/FLTConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${FLT_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(FLTConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FLTConfig.cmake" @ONLY)
# ... for both
configure_file(FLTConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/FLTConfigVersion.cmake" @ONLY)

# Install the FltConfig.cmake and FltConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FLTConfig.cmake"
  "${PROJECT_BINARY_DIR}/FLTConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

# Install the export set for use with the install-tree
install(EXPORT FLTTargets DESTINATION
  "${INSTALL_CMAKE_DIR}")
