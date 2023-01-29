add_executable(main main.cpp)

add_subdirectory(graphutils)

list(APPEND MAIN_LIBS Graphutils)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/graphutils") 
target_link_libraries(main PUBLIC ${MAIN_LIBS})

add_subdirectory(tests)
