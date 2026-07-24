#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <utility>
#include <memory>
#include <string>
#include <ortools/linear_solver/linear_solver.h>

namespace aal {

class Simplex {
public:
    int m, n;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> obj_row;

    Simplex(const std::vector<std::vector<double>>& A_in, const std::vector<double>& b_in, const std::vector<double>& c_in) {
        m = A_in.size();
        n = m > 0 ? A_in[0].size() : 0;
        A = A_in;
        b = b_in;
        c = c_in;
    }

    std::pair<std::vector<double>, double> solve() {
        std::unique_ptr<operations_research::MPSolver> solver(
            operations_research::MPSolver::CreateSolver("GLOP"));
        if (!solver) {
            std::cerr << "GLOP solver not available.\n";
            return {{}, std::numeric_limits<double>::infinity()};
        }

        const double infinity = solver->infinity();

        std::vector<operations_research::MPVariable*> y_vars(n);
        for (int j = 0; j < n; ++j) {
            y_vars[j] = solver->MakeNumVar(0.0, infinity, "y_" + std::to_string(j));
        }

        std::vector<operations_research::MPConstraint*> constraints(m);
        for (int i = 0; i < m; ++i) {
            constraints[i] = solver->MakeRowConstraint(-infinity, b[i], "c_" + std::to_string(i));
            for (int j = 0; j < n; ++j) {
                constraints[i]->SetCoefficient(y_vars[j], A[i][j]);
            }
        }

        operations_research::MPObjective* const objective = solver->MutableObjective();
        for (int j = 0; j < n; ++j) {
            objective->SetCoefficient(y_vars[j], c[j]);
        }
        objective->SetMaximization();

        const operations_research::MPSolver::ResultStatus result_status = solver->Solve();

        if (result_status != operations_research::MPSolver::OPTIMAL) {
            return {{}, std::numeric_limits<double>::infinity()};
        }

        std::vector<double> y_sol(n);
        for (int j = 0; j < n; ++j) {
            y_sol[j] = y_vars[j]->solution_value();
        }

        obj_row.assign(n + m + 1, 0.0);
        for (int j = 0; j < n; ++j) {
            obj_row[j] = y_vars[j]->reduced_cost();
        }
        for (int i = 0; i < m; ++i) {
            obj_row[n + i] = constraints[i]->dual_value();
        }
        obj_row.back() = objective->Value();

        return {y_sol, objective->Value()};
    }
};

} // namespace aal
