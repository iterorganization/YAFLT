# RKF45 dependency
set(RKF45_SOURCES
    rkf45.cpp
)

# Bicubic dependency
set(BICUBIC_SOURCES
    bicubic.cpp
)

# FLT main library
set(FLT_SOURCES
    flt.cpp
    accell_embree.cpp
)

# All the header files, which we install in the include directory of the
# installation directory
set(FLT_HEADERS
    rkf45.hpp
    flt.hpp
    accell_embree.hpp
    bicubic.hpp
    )

# Create the static library and shared library
add_library(flt_static STATIC 
    ${BICUBIC_SOURCES} ${RKF45_SOURCES} ${FLT_SOURCES})
add_library(flt SHARED 
    ${BICUBIC_SOURCES} ${RKF45_SOURCES} ${FLT_SOURCES})

# Link embree and add include directories for dependencies.
foreach(target flt_static flt)
    target_link_libraries(${target} PUBLIC ${EMBREE_LIBRARY})
    target_include_directories(${target} PUBLIC ${EMBREE_INCLUDE_DIRS})
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
endforeach(target)

if(NOT MSVC)
    # Add the headers to the PUBLIC_HEADER property so that they get installed.
    # set_target_properties(flt_static PROPERTIES PUBLIC_HEADER "${FLT_HEADERS}")
    # Rename the static library so that in the end we will have a flt.a and flt.so
    # and not flt_static.a and flt.so.
    set_target_properties(flt_static PROPERTIES OUTPUT_NAME "flt")
    set_target_properties(flt PROPERTIES CXX_VISIBILITY_PRESET hidden)
else()
    target_compile_definitions(flt PRIVATE "flt_EXPORTS")
    # PUBLIC required when linking the static library.
    target_compile_definitions(flt_static PUBLIC "flt_STATIC")
endif()
# Set the list of header files to be installed in the end.
set_target_properties(flt_static PROPERTIES PUBLIC_HEADER "${FLT_HEADERS}")

# Mark libraries for installation. In case of flt we also install the public
# headers.
install(TARGETS flt 
    EXPORT FLTTargets
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
    RUNTIME DESTINATION "${INSTALL_BIN_DIR}"
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}" 
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")

install(TARGETS flt_static
    EXPORT FLTTargets
    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
    PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}")
