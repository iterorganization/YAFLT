#ALGLIB sources
set(ALGLIB_SOURCES
    ../include/alglib/alglibinternal.cpp
    ../include/alglib/alglibmisc.cpp
    ../include/alglib/ap.cpp
    ../include/alglib/dataanalysis.cpp
    ../include/alglib/diffequations.cpp
    ../include/alglib/fasttransforms.cpp
    ../include/alglib/integration.cpp
    ../include/alglib/interpolation.cpp
    ../include/alglib/linalg.cpp
    ../include/alglib/optimization.cpp
    ../include/alglib/solvers.cpp
    ../include/alglib/specialfunctions.cpp
    ../include/alglib/statistics.cpp
)

add_library(alglib_static STATIC ${ALGLIB_SOURCES})

set(BICUBIC_SOURCES
    ../../src/bicubic.cpp
)

add_library(bicubic_scalable STATIC ${BICUBIC_SOURCES})

# add the executable
foreach(target 1_accuracy 3_performance 4_derivative_eval 6_check_values 7_test_rkf45_stability 8_simple_test_case)
    add_executable(${target} ${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp)
    target_link_libraries(${target} flt_static)

    # For linking example header data.
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/alglib)
    target_link_libraries(${target} alglib_static)
    # install(TARGETS ${target} RUNTIME DESTINATION bin)
endforeach()


# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    add_executable(9_omp_performance ${CMAKE_CURRENT_SOURCE_DIR}/9_omp_performance.cpp)
    target_link_libraries(9_omp_performance PUBLIC flt_static)
    # For linking example header data.
    target_include_directories(9_omp_performance PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_link_libraries(9_omp_performance PUBLIC bicubic_scalable)
    target_link_libraries(9_omp_performance PUBLIC OpenMP::OpenMP_CXX)

endif(OPENMP_FOUND)