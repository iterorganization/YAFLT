# Bicubic tests
add_executable(test_1_accuracy ${CMAKE_CURRENT_SOURCE_DIR}/test_1_accuracy.cpp)
target_link_libraries(test_1_accuracy flt_static)

add_executable(test_2_serial_performance ${CMAKE_CURRENT_SOURCE_DIR}/test_2_serial_performance.cpp)
target_link_libraries(test_2_serial_performance flt_static)
target_include_directories(test_2_serial_performance PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    add_executable(test_3_omp_performance ${CMAKE_CURRENT_SOURCE_DIR}/test_3_omp_performance.cpp)
    target_link_libraries(test_3_omp_performance flt_static)
    target_include_directories(test_3_omp_performance PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
endif(OPENMP_FOUND)

add_executable(test_4_finite_difference_check ${CMAKE_CURRENT_SOURCE_DIR}/test_4_finite_difference_check.cpp)
target_link_libraries(test_4_finite_difference_check flt_static)
target_include_directories(test_4_finite_difference_check PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)
