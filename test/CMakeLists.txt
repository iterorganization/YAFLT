set(INTERPOLATION_TESTS
    interpolation/accuracy.cpp
    interpolation/serial_performance.cpp
    interpolation/omp_performance.cpp
    interpolation/finite_difference_check.cpp
)

set(OBJECT_TESTS
    objects/create_and_delete_objects.cpp
)

set(RKF45_TESTS
    rkf45/test_stability.cpp
    rkf45/test_rkf45_method.cpp
    rkf45/test_rkf45_method_stability.cpp
)

set(EMBREE_TESTS
    embree/load_mesh.cpp
    embree/test_intersection.cpp
)

foreach(target_src ${INTERPOLATION_TESTS} ${OBJECT_TESTS} ${RKF45_TESTS} ${EMBREE_TESTS})
    get_filename_component(category ${target_src} DIRECTORY)
    get_filename_component(target_name ${target_src} NAME_WE) 
    set(executable ${category}_${target_name})
    message("-- Creating test ${executable} ${target_src}")

    add_executable(${executable} ${target_src})
    target_link_libraries(${executable} flt_static)
    target_include_directories(${executable} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(${executable} PRIVATE $<TARGET_PROPERTY:flt_static,INCLUDE_DIRECTORIES>)
    add_test(NAME "${category}/${target_name}" COMMAND $<TARGET_FILE:${executable}>)
endforeach(target_src)