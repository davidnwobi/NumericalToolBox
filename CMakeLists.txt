cmake_minimum_required(VERSION 3.12)
project(NumLibUP)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(BOOST_ROOT "../Packages/boost_1_83_0")
add_subdirectory(googletest)

find_package(Boost REQUIRED)


add_library(NumLib
    src/gnuplot-iostream.h
    src/linodesolver.h
    src/LSF.h
    src/matrixtemplate.h
    src/matrixtoolbox.h
    src/numerical_functions.h
    src/print.h
    src/vectortemplate.h
)
set_target_properties(NumLib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(NumLib PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(NumLib PUBLIC ${Boost_LIBRARIES})

add_executable(main 
    main.cpp
)


target_link_libraries(main PRIVATE NumLib)
add_compile_options(-O3 -funroll-loops)
enable_testing()
add_subdirectory(tests)