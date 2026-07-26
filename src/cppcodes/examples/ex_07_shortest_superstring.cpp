/**
 * Chapter 7: Shortest Superstring (Greedy Cycle Cover)
 * 
 * Theory:
 *   Given a set of strings, find the shortest string containing all of them as substrings.
 *   This is NP-hard. We implement the greedy overlap heuristic and the 4-approx / 3-approx
 *   cycle cover based algorithms.
 */
#include "../approx_algorithms.hpp"
#include <iostream>

int main() {
    using namespace aal;
    std::cout << "=== Chapter 7: Shortest Superstring ===\n";

    std::vector<std::string> strings = {"CAT", "ATC", "TCA"};
    std::cout << "  Input Strings: [\"CAT\", \"ATC\", \"TCA\"]\n";

    std::string greedy = greedy_superstring(strings);
    std::string s_4ap = shortest_superstring_4approx(strings);
    std::string s_3ap = shortest_superstring_3approx(strings);

    std::cout << "  Greedy Heuristic Superstring:      " << greedy << " (len " << greedy.length() << ")\n";
    std::cout << "  Cycle-Cover 4-approx Superstring:   " << s_4ap << " (len " << s_4ap.length() << ")\n";
    std::cout << "  Cycle-Cover 3-approx Superstring:   " << s_3ap << " (len " << s_3ap.length() << ")\n\n";
    return 0;
}
