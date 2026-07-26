#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <map>
#include <iomanip>

namespace aal {

double compute_tour_cost(const std::vector<int>& tour, const std::vector<std::pair<double, double>>& points) {
    double cost = 0.0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        int u = tour[i];
        int v = tour[i + 1];
        cost += std::hypot(points[u].first - points[v].first, points[u].second - points[v].second);
    }
    return cost;
}

std::pair<std::vector<int>, double> held_karp_tsp(const std::vector<std::pair<double, double>>& points) {
    int n = points.size();
    if (n == 0) return {{}, 0.0};
    if (n == 1) return {{0}, 0.0};
    
    std::vector<std::vector<double>> dist(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dist[i][j] = std::hypot(points[i].first - points[j].first, points[i].second - points[j].second);
        }
    }
    
    std::map<std::pair<int, int>, double> memo;
    std::map<std::pair<int, int>, int> parent;
    
    auto tsp_solve = [&](auto& self, int mask, int last) -> double {
        if (mask == (1 << n) - 1) {
            return dist[last][0];
        }
        std::pair<int, int> state = {mask, last};
        if (memo.count(state)) return memo[state];
        
        double best = 1e9;
        int best_next = -1;
        for (int nxt = 0; nxt < n; ++nxt) {
            if (!(mask & (1 << nxt))) {
                double cost = dist[last][nxt] + self(self, mask | (1 << nxt), nxt);
                if (cost < best) {
                    best = cost;
                    best_next = nxt;
                }
            }
        }
        memo[state] = best;
        parent[state] = best_next;
        return best;
    };
    
    double opt_cost = tsp_solve(tsp_solve, 1, 0);
    
    std::vector<int> path = {0};
    int mask = 1;
    int curr = 0;
    for (int i = 0; i < n - 1; ++i) {
        int nxt = parent[{mask, curr}];
        path.push_back(nxt);
        mask |= (1 << nxt);
        curr = nxt;
    }
    path.push_back(0);
    
    return {path, opt_cost};
}

std::vector<int> merge_tours(const std::vector<int>& tour1, const std::vector<int>& tour2, const std::vector<std::pair<double, double>>& points) {
    if (tour1.empty()) return tour2;
    if (tour2.empty()) return tour1;
    
    auto dist = [&](int u, int v) {
        return std::hypot(points[u].first - points[v].first, points[u].second - points[v].second);
    };
    
    if (tour2.size() == 1 || (tour2.size() == 2 && tour2[0] == tour2[1])) {
        int u = tour2[0];
        double best_cost_diff = 1e9;
        int best_idx = -1;
        for (size_t i = 0; i < tour1.size() - 1; ++i) {
            int u1 = tour1[i];
            int v1 = tour1[i + 1];
            double diff = dist(u1, u) + dist(u, v1) - dist(u1, v1);
            if (diff < best_cost_diff) {
                best_cost_diff = diff;
                best_idx = i;
            }
        }
        std::vector<int> merged;
        for (int i = 0; i <= best_idx; ++i) merged.push_back(tour1[i]);
        merged.push_back(u);
        for (size_t i = best_idx + 1; i < tour1.size(); ++i) merged.push_back(tour1[i]);
        return merged;
    }
    
    if (tour1.size() == 1 || (tour1.size() == 2 && tour1[0] == tour1[1])) {
        return merge_tours(tour2, tour1, points);
    }
    
    std::vector<int> t1 = tour1;
    if (t1.back() != t1.front()) t1.push_back(t1.front());
    std::vector<int> t2 = tour2;
    if (t2.back() != t2.front()) t2.push_back(t2.front());
    
    int n1 = t1.size() - 1;
    int n2 = t2.size() - 1;
    double best_cost_diff = 1e9;
    int best_i = -1, best_j = -1;
    bool best_reverse = false;
    
    for (int i = 0; i < n1; ++i) {
        int u1 = t1[i];
        int v1 = t1[i + 1];
        for (int j = 0; j < n2; ++j) {
            int u2 = t2[j];
            int v2 = t2[j + 1];
            
            double diff1 = dist(u1, u2) + dist(v1, v2) - dist(u1, v1) - dist(u2, v2);
            double diff2 = dist(u1, v2) + dist(v1, u2) - dist(u1, v1) - dist(u2, v2);
            
            if (diff1 < best_cost_diff) {
                best_cost_diff = diff1;
                best_i = i; best_j = j; best_reverse = false;
            }
            if (diff2 < best_cost_diff) {
                best_cost_diff = diff2;
                best_i = i; best_j = j; best_reverse = true;
            }
        }
    }
    
    if (best_i == -1) {
        std::vector<int> res = t1;
        res.insert(res.end(), t2.begin() + 1, t2.end());
        return res;
    }
    
    std::vector<int> part1(t1.begin(), t1.begin() + best_i + 1);
    std::vector<int> t2_rotated;
    t2_rotated.insert(t2_rotated.end(), t2.begin() + best_j, t2.end() - 1);
    t2_rotated.insert(t2_rotated.end(), t2.begin(), t2.begin() + best_j);
    
    if (best_reverse) {
        std::reverse(t2_rotated.begin(), t2_rotated.end());
    }
    
    std::vector<int> merged = part1;
    merged.insert(merged.end(), t2_rotated.begin(), t2_rotated.end());
    merged.insert(merged.end(), t1.begin() + best_i + 1, t1.end());
    
    if (merged.back() != merged.front()) merged.push_back(merged.front());
    return merged;
}

std::vector<int> quadtree_tsp_solve(const std::vector<int>& point_indices, const std::vector<std::pair<double, double>>& points) {
    int n = point_indices.size();
    if (n == 0) return {};
    if (n <= 3) {
        std::vector<std::pair<double, double>> local_points;
        for (int idx : point_indices) local_points.push_back(points[idx]);
        auto res = held_karp_tsp(local_points);
        std::vector<int> global_tour;
        for (int idx : res.first) global_tour.push_back(point_indices[idx]);
        return global_tour;
    }
    
    double min_x = 1e9, max_x = -1e9, min_y = 1e9, max_y = -1e9;
    for (int idx : point_indices) {
        min_x = std::min(min_x, points[idx].first);
        max_x = std::max(max_x, points[idx].first);
        min_y = std::min(min_y, points[idx].second);
        max_y = std::max(max_y, points[idx].second);
    }
    
    double mid_x = (min_x + max_x) / 2;
    double mid_y = (min_y + max_y) / 2;
    
    if (max_x - min_x < 1e-9 && max_y - min_y < 1e-9) {
        std::vector<int> tour = point_indices;
        tour.push_back(point_indices[0]);
        return tour;
    }
    
    std::vector<int> q1, q2, q3, q4;
    for (int idx : point_indices) {
        double x = points[idx].first, y = points[idx].second;
        if (x <= mid_x) {
            if (y <= mid_y) q1.push_back(idx);
            else q2.push_back(idx);
        } else {
            if (y <= mid_y) q3.push_back(idx);
            else q4.push_back(idx);
        }
    }
    
    if (q1.size() == n || q2.size() == n || q3.size() == n || q4.size() == n) {
        int half = n / 2;
        q1 = std::vector<int>(point_indices.begin(), point_indices.begin() + half);
        q2 = std::vector<int>(point_indices.begin() + half, point_indices.end());
        q3.clear(); q4.clear();
    }
    
    auto t1 = quadtree_tsp_solve(q1, points);
    auto t2 = quadtree_tsp_solve(q2, points);
    auto t3 = quadtree_tsp_solve(q3, points);
    auto t4 = quadtree_tsp_solve(q4, points);
    
    std::vector<int> merged = t1;
    if (!t2.empty()) merged = merge_tours(merged, t2, points);
    if (!t3.empty()) merged = merge_tours(merged, t3, points);
    if (!t4.empty()) merged = merge_tours(merged, t4, points);
    
    return merged;
}

std::pair<std::vector<int>, double> quadtree_tsp(const std::vector<std::pair<double, double>>& points) {
    std::vector<int> indices(points.size());
    for (size_t i = 0; i < points.size(); ++i) indices[i] = i;
    
    auto tour = quadtree_tsp_solve(indices, points);
    double cost = compute_tour_cost(tour, points);
    return {tour, cost};
}

} // namespace aal

using namespace aal;

void demo_euclidean_tsp() {
    std::cout << "============================================================\n";
    std::cout << "Chapter 11: Euclidean TSP Heuristics\n";
    std::cout << "============================================================\n";
    
    std::vector<std::pair<double, double>> points = {
        {0.0, 0.0}, {1.0, 4.0}, {3.0, 1.0}, {4.0, 3.0},
        {1.0, 1.0}, {3.0, 4.0}, {5.0, 0.0}, {5.5, 4.5}
    };
    
    std::cout << "\n1. Input 2D Points (n=" << points.size() << "):\n";
    for (size_t idx = 0; idx < points.size(); ++idx) {
        std::cout << "  Point " << idx << ": (" << points[idx].first << ", " << points[idx].second << ")\n";
    }
    
    auto [opt_tour, opt_cost] = held_karp_tsp(points);
    std::cout << "\nExact Held-Karp Tour: [";
    for (size_t i = 0; i < opt_tour.size(); ++i) std::cout << opt_tour[i] << (i + 1 == opt_tour.size() ? "" : ", ");
    std::cout << "]\nExact Optimal Cost:    " << std::fixed << std::setprecision(4) << opt_cost << "\n";
    
    auto [qt_tour, qt_cost] = quadtree_tsp(points);
    std::cout << "\nQuadtree Heuristic Tour: [";
    for (size_t i = 0; i < qt_tour.size(); ++i) std::cout << qt_tour[i] << (i + 1 == qt_tour.size() ? "" : ", ");
    std::cout << "]\nHeuristic Tour Cost:     " << std::setprecision(4) << qt_cost << "\n";
    std::cout << "Approximation Ratio:     " << (qt_cost / opt_cost) << " (Theoretical: 1 + eps)\n";
}
