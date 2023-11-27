include(GoogleTest)

set(TEST_FILES 
    CVectorTest.cpp
    MatrixTest.cpp
    MatrixtoolboxTest.cpp
)

add_executable(NumLibTests ${TEST_FILES})

target_link_libraries(NumLibTests PRIVATE gtest_main NumLib)

target_include_directories(NumLibTests PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
    ${Boost_INCLUDE_DIRS}
)

gtest_discover_tests(NumLibTests)
