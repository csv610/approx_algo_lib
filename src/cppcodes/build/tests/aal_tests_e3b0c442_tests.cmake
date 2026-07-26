add_test([=[VertexCoverTest.EmptyGraph]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=VertexCoverTest.EmptyGraph]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[VertexCoverTest.EmptyGraph]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_vertex_cover.cpp:4]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
add_test([=[VertexCoverTest.Path5approx]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=VertexCoverTest.Path5approx]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[VertexCoverTest.Path5approx]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_vertex_cover.cpp:11]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
add_test([=[SetCoverTest.BasicCover]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=SetCoverTest.BasicCover]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[SetCoverTest.BasicCover]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_set_cover.cpp:4]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
add_test([=[KnapsackTest.ExactAndFPTAS]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=KnapsackTest.ExactAndFPTAS]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[KnapsackTest.ExactAndFPTAS]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_knapsack.cpp:4]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
add_test([=[SimplexTest.SimpleLP]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=SimplexTest.SimpleLP]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[SimplexTest.SimpleLP]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_simplex.cpp:4]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
add_test([=[SteinerTest.MetricTSP]=]  /Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests/aal_tests [==[--gtest_filter=SteinerTest.MetricTSP]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[SteinerTest.MetricTSP]=]
  PROPERTIES
    
    DEF_SOURCE_LINE [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/tests/test_steiner.cpp:4]==]
    WORKING_DIRECTORY [==[/Users/csv610/Projects/MyBooks/ApproxAlgo/approx_algorithms_cpp/build/tests]==]
    SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==]
    
)
set(aal_tests_TESTS [==[VertexCoverTest.EmptyGraph]==] [==[VertexCoverTest.Path5approx]==] [==[SetCoverTest.BasicCover]==] [==[KnapsackTest.ExactAndFPTAS]==] [==[SimplexTest.SimpleLP]==] [==[SteinerTest.MetricTSP]==])
