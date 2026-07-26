/**
 * Chapter 10: Makespan PTAS
 * 
 * Theory:
 *   Scheduling n independent jobs on m machines to minimize makespan (P || C_max).
 *   The PTAS schedules large jobs exactly via DP and schedules small jobs
 *   greedily (List Scheduling), guaranteeing makespan <= (1 + eps)*OPT.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 10: Makespan PTAS ===\n";

    std::vector<double> jobs = {2.0, 3.0, 4.0, 5.0, 6.0};
    int m = 3;
    double eps = 0.25;

    std::cout << "  Jobs: [2.0, 3.0, 4.0, 5.0, 6.0], Machines: " << m << ", Epsilon: " << eps << "\n";

    auto schedule = makespan_ptas(jobs, m, eps);

    std::cout << "  PTAS Scheduled Machines Count: " << schedule.size() << "\n\n";
    return 0;
}
