# Alglib dependency
set(ALGLIB_SOURCES
    alglib/alglibinternal.cpp
    alglib/alglibmisc.cpp
    alglib/ap.cpp
    alglib/dataanalysis.cpp
    alglib/diffequations.cpp
    alglib/fasttransforms.cpp
    alglib/integration.cpp
    alglib/interpolation.cpp
    alglib/linalg.cpp
    alglib/optimization.cpp
    alglib/solvers.cpp
    alglib/specialfunctions.cpp
    alglib/statistics.cpp
)

# RKF45 dependency
set(RKF45_SOURCES
    rkf45.cpp
)

# FLT main library
set(FLT_SOURCES
    flt.cpp
    accell_embree.cpp
)

# All the header files, which we install in the include directory of the
# installation directory
set(FLT_HEADERS
    alglib/integration.h
    alglib/optimization.h
    alglib/stdafx.h
    alglib/diffequations.h
    alglib/specialfunctions.h
    alglib/alglibinternal.h
    alglib/alglibmisc.h
    alglib/statistics.h
    alglib/linalg.h
    alglib/dataanalysis.h
    alglib/fasttransforms.h
    alglib/solvers.h
    alglib/interpolation.h
    alglib/ap.h
    rkf45.hpp
    flt.hpp
    accell_embree.hpp
    )

# Create the static library and shared library
add_library(flt_static STATIC ${ALGLIB_SOURCES} ${RKF45_SOURCES} ${FLT_SOURCES})
add_library(flt SHARED ${ALGLIB_SOURCES} ${RKF45_SOURCES} ${FLT_SOURCES})

# Link embree and add include directories for dependencies.
foreach(target flt_static flt)
    target_link_libraries(${target} ${EMBREE_LIBRARY})
    target_include_directories(${target} PUBLIC ${EMBREE_INCLUDE_DIRS})
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/alglib)
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/rkf45)
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endforeach(target)

# Add the headers to the PUBLIC_HEADER property so that they get installed.
# set_target_properties(flt_static PROPERTIES PUBLIC_HEADER "${FLT_HEADERS}")
# Rename the static library so that in the end we will have a flt.a and flt.so
# and not flt_static.a and flt.so.
set_target_properties(flt_static PROPERTIES OUTPUT_NAME "flt")
# Set the list of header files to be installed in the end.
set_target_properties(flt_static PROPERTIES PUBLIC_HEADER "${FLT_HEADERS}")

# Mark libraries for installation. In case of flt we also install the public
# headers.
install(TARGETS flt LIBRARY DESTINATION lib PUBLIC_HEADER DESTINATION include/flt RUNTIME DESTINATION bin)
install(TARGETS flt_static DESTINATION lib PUBLIC_HEADER DESTINATION include/flt)
