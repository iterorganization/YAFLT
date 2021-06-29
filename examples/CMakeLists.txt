# add the executable
foreach(target flt1 flt2 flt3 embreeTest1)
    add_executable(${target} ${PROJECT_SOURCE_DIR}/examples/${target}.cpp)
    target_link_libraries(${target} flt_static)

    # For linking example header data.
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

    install(TARGETS ${target} RUNTIME DESTINATION bin)
endforeach()


