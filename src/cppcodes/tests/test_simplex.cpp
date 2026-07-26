#include <gtest/gtest.h>
#include "../simplex.hpp"

TEST(SimplexTest, SimpleLP) {
    using namespace aal;
    
    // Maximize 3x_1 + 4x_2
    // Subject to:
    //   x_1 + 2x_2 <= 8
    //   3x_1 + 2x_2 <= 12
    //   x_1, x_2 >= 0
    std::vector<std::vector<double>> A = {
        {1.0, 2.0},
        {3.0, 2.0}
    };
    std::vector<double> b = {8.0, 12.0};
    std::vector<double> c = {3.0, 4.0};

    Simplex solver(A, b, c);
    auto [x, opt] = solver.solve();

    EXPECT_NEAR(opt, 19.0, 1e-6); // Max objective value is 19.0 at (2, 3)
    ASSERT_EQ(x.size(), 2);
    EXPECT_NEAR(x[0], 2.0, 1e-6);
    EXPECT_NEAR(x[1], 3.0, 1e-6);
}
