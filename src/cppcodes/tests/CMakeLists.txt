cmake_minimum_required(VERSION 3.20)

find_package(GTest REQUIRED)

add_executable(aal_tests
    test_vertex_cover.cpp
    test_set_cover.cpp
    test_knapsack.cpp
    test_simplex.cpp
    test_steiner.cpp
)

target_link_libraries(aal_tests PRIVATE approx_algo GTest::gtest GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(aal_tests)
