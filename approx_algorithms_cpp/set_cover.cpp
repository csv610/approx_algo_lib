#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cmath>

namespace aal {

using Universe = std::set<int>;

std::pair<std::vector<int>, double> greedy_set_cover(
    const Universe& universe, 
    const std::map<int, std::set<int>>& sets, 
    const std::map<int, double>& costs) {
    
    std::set<int> covered;
    std::vector<int> picked_sets;
    double total_cost = 0.0;
    std::map<int, double> prices;
    
    while (covered != universe) {
        int best_set = -1;
        double best_cost_effectiveness = std::numeric_limits<double>::infinity();
        
        for (const auto& [sid, s] : sets) {
            std::set<int> new_elements;
            std::set_difference(s.begin(), s.end(), covered.begin(), covered.end(),
                                std::inserter(new_elements, new_elements.begin()));
            
            if (new_elements.empty()) continue;
            
            double cost_effectiveness = costs.at(sid) / new_elements.size();
            if (cost_effectiveness < best_cost_effectiveness) {
                best_cost_effectiveness = cost_effectiveness;
                best_set = sid;
            }
        }
        
        if (best_set == -1) break;
        
        std::set<int> new_elements;
        const auto& s = sets.at(best_set);
        std::set_difference(s.begin(), s.end(), covered.begin(), covered.end(),
                            std::inserter(new_elements, new_elements.begin()));
        
        for (int e : new_elements) {
            prices[e] = best_cost_effectiveness;
        }
        
        covered.insert(new_elements.begin(), new_elements.end());
        picked_sets.push_back(best_set);
        total_cost += costs.at(best_set);
    }
    
    return {picked_sets, total_cost};
}

double harmonic_number(int n) {
    double sum = 0.0;
    for (int i = 1; i <= n; ++i) sum += 1.0 / i;
    return sum;
}

struct SetCoverInstance {
    Universe universe;
    std::map<int, std::set<int>> sets;
    std::map<int, double> costs;
};

SetCoverInstance set_cover_tight_example(int n) {
    SetCoverInstance instance;
    for (int i = 1; i <= n; ++i) {
        instance.universe.insert(i);
        instance.sets[0].insert(i);
    }
    instance.costs[0] = 1.0 + 1e-9;
    
    for (int i = 1; i <= n; ++i) {
        instance.sets[i] = {i};
        instance.costs[i] = 1.0 / i;
    }
    return instance;
}

std::pair<std::vector<int>, double> set_cover_exact_bruteforce(
    const Universe& universe, 
    const std::map<int, std::set<int>>& sets, 
    const std::map<int, double>& costs) {
    
    std::vector<int> set_ids;
    for (const auto& kv : sets) set_ids.push_back(kv.first);
    
    std::vector<int> best_sets = set_ids;
    double best_cost = std::numeric_limits<double>::infinity();
    
    int num_sets = set_ids.size();
    for (int mask = 0; mask < (1 << num_sets); ++mask) {
        std::set<int> cover;
        double cost = 0.0;
        for (int i = 0; i < num_sets; ++i) {
            if (mask & (1 << i)) {
                int sid = set_ids[i];
                cover.insert(sets.at(sid).begin(), sets.at(sid).end());
                cost += costs.at(sid);
            }
        }
        if (cover == universe && cost < best_cost) {
            best_cost = cost;
            best_sets.clear();
            for (int i = 0; i < num_sets; ++i) {
                if (mask & (1 << i)) {
                    best_sets.push_back(set_ids[i]);
                }
            }
        }
    }
    return {best_sets, best_cost};
}

} // namespace aal

void demo_set_cover() {
    using namespace aal;
    std::cout << "============================================================\n";
    std::cout << "Chapter 2: Set Cover - Greedy Algorithm (H_n approximation)\n";
    std::cout << "============================================================\n";
    
    std::cout << "\n1. Tight Example (Vazirani Example 2.5)\n";
    for (int n : {5, 10, 20}) {
        auto instance = set_cover_tight_example(n);
        auto [greedy_sets, greedy_cost] = greedy_set_cover(instance.universe, instance.sets, instance.costs);
        auto [optimal_sets, optimal_cost] = set_cover_exact_bruteforce(instance.universe, instance.sets, instance.costs);
        double Hn = harmonic_number(n);
        std::cout << "  n=" << n << ": Greedy=" << greedy_cost << ", Opt=" << optimal_cost 
                  << ", Ratio=" << greedy_cost / optimal_cost << ", H_n=" << Hn << "\n";
    }
    
    std::cout << "\n2. Practical Example: Feature Selection\n";
    Universe universe;
    for (int i = 0; i < 20; ++i) universe.insert(i);
    std::map<int, std::set<int>> sets = {
        {0, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}},
        {1, {10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {2, {0, 2, 4, 6, 8, 10, 12, 14, 16, 18}},
        {3, {1, 3, 5, 7, 9, 11, 13, 15, 17, 19}},
        {4, {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}}
    };
    std::map<int, double> costs = {{0, 10.0}, {1, 10.0}, {2, 6.0}, {3, 6.0}, {4, 8.0}};
    
    auto [greedy_sets, greedy_cost] = greedy_set_cover(universe, sets, costs);
    std::cout << "  Universe size: " << universe.size() << "\n";
    std::cout << "  Available sets: " << sets.size() << "\n";
    std::cout << "  Greedy picked: [";
    for (int s : greedy_sets) std::cout << s << " ";
    std::cout << "], cost=" << greedy_cost << "\n";
}
