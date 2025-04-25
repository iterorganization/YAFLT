# add the executable
foreach(target flt1 flt2 flt3 inres1_example parallel_speed_testing)
    add_executable(${target} ${PROJECT_SOURCE_DIR}/examples/${target}.cpp)
    target_link_libraries(${target} flt_static)

    # For linking example header data.
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${target} PRIVATE $<TARGET_PROPERTY:flt_static,INCLUDE_DIRECTORIES>)
    if(MSVC)
        # Increase the stack size as the data headers cause a segfault without
        # increased stack.
        set_target_properties(${target} PROPERTIES LINK_FLAGS "/STACK:8388608")
    endif()
endforeach()
