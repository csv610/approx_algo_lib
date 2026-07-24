cmake_minimum_required(VERSION 3.20)
project(ApproxAlgorithmsCPP CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add all warning flags
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -O3)
endif()

find_package(ortools REQUIRED)

add_library(approx_algo STATIC
    ch01_intro.cpp
    ch02_set_cover.cpp
    ch03_steiner_tsp.cpp
    ch04_multiway_kcut.cpp
    ch05_kcenter.cpp
    ch06_feedback_vertex_set.cpp
    ch07_shortest_superstring.cpp
    ch08_knapsack.cpp
    ch09_bin_packing.cpp
    ch10_makespan.cpp
    ch11_euclidean_tsp.cpp
    ch12_14_lp_algorithms.cpp
    ch13_knapsack.cpp
    ch15_weighted_vertex_cover_pd.cpp
    ch16_randomized_rounding.cpp
    ch17_chernoff_bounds.cpp
    ch18_sdp_maxcut.cpp
    ch19_multiway_rounding.cpp
    ch21_steiner_forest.cpp
    ch22_steiner_network.cpp
    ch23_primal_dual_fvs.cpp
    ch24_facility_location.cpp
    ch26_sdp_sat.cpp
    ch30_tree_multicut.cpp
)
target_link_libraries(approx_algo PUBLIC ortools::ortools)

add_executable(approx_algo_demos main.cpp)
target_link_libraries(approx_algo_demos PRIVATE approx_algo)

add_subdirectory(examples)
