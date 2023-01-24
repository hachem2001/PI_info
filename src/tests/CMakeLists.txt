add_executable(test1 test1.cpp)

#add_subdirectory(../graphutils)
#list(APPEND MAIN_LIBS Graphutils)

target_link_libraries(test1 PRIVATE ${MAIN_LIBS})
