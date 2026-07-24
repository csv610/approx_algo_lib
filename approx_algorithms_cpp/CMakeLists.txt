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
    intro.cpp
    set_cover.cpp
    steiner_tsp.cpp
    multiway_kcut.cpp
    kcenter.cpp
    feedback_vertex_set.cpp
    shortest_superstring.cpp
    knapsack.cpp
    bin_packing.cpp
    makespan.cpp
    euclidean_tsp.cpp
    lp_algorithms.cpp
    knapsack_ch13.cpp
    weighted_vertex_cover_pd.cpp
    randomized_rounding.cpp
    chernoff_bounds.cpp
    sdp_maxcut.cpp
    multiway_rounding.cpp
    steiner_forest.cpp
    steiner_network.cpp
    primal_dual_fvs.cpp
    facility_location.cpp
    sdp_sat.cpp
    tree_multicut.cpp
)
target_link_libraries(approx_algo PUBLIC ortools::ortools)

add_executable(approx_algo_demos main.cpp)
target_link_libraries(approx_algo_demos PRIVATE approx_algo)

add_subdirectory(examples)
