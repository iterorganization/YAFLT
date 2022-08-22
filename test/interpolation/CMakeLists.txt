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

# add the executable
foreach(target 1_accuracy 2_accuracy_with_flt 3_performance 4_derivative_eval 5_derivative_eval_with_flt 6_check_values 7_test_rkf45_stability)
    add_executable(${target} ${CMAKE_CURRENT_SOURCE_DIR}/${target}.cpp)
    target_link_libraries(${target} flt_static)

    # For linking example header data.
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/alglib)
    target_link_libraries(${target} alglib_static)
    # install(TARGETS ${target} RUNTIME DESTINATION bin)
endforeach()

