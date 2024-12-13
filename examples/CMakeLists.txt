# add the executable
foreach(target flt1 flt2 flt3 inres1_example)
    add_executable(${target} ${PROJECT_SOURCE_DIR}/examples/${target}.cpp)
    target_link_libraries(${target} flt_static)

    # For linking example header data.
    target_include_directories(${target} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_include_directories(${target} PRIVATE $<TARGET_PROPERTY:flt_static,INCLUDE_DIRECTORIES>)

    install(TARGETS ${target} RUNTIME DESTINATION ${INSTALL_BIN_DIR})
endforeach()


