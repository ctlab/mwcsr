//
// Created by Alexander Loboda on 9/15/21.
//

#include "include/primal_heuristic.h"

#include <functional>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <limits>

namespace {
class Update {
    size_t e;
    size_t iteration;
    double distance;
public:
    Update(size_t e, size_t iteration, double d) :e(e), iteration(iteration), distance(d) {}
    size_t edge() const {
        return e;
    }

    size_t it() const {
        return iteration;
    }

    double d() const {
        return distance;
    }

    bool operator<(const Update& upd) const {
        return distance > upd.distance;
    }
};

}

namespace relax {

PrimalHeuristic::PrimalHeuristic(mwcsr::Graph g, std::function<double(size_t)> wf, std::vector<size_t> active,
                                 std::vector<bool> current)
    :g(g), weight(wf), active(active), current(current) {
    for(auto e: active) {
        shortest_distance[e] = std::numeric_limits<double>::infinity();
    }
}

Solution PrimalHeuristic::run_heuristic() {
    std::sort(active.begin(), active.end(), [this](size_t l, size_t r) {
        if (current.at(l) != current.at(r)) {
            return current.at(l) > current.at(r);
        }
        return wei(l) > wei(r);
    });
    for (size_t e: active) {
        if (visit.find(e) == visit.end()) {
            run_from_point(e);
        }
    }
    return best;
}

void PrimalHeuristic::run_from_point(size_t e) {
    std::unordered_set<size_t> signals;
    std::unordered_set<size_t> in_solution;
    std::priority_queue<Update> q;
    shortest_distance[e] = 0.0;
    previous[e] = -1;
    q.emplace(e, iteration, 0.0);

    Solution solution;
    while (!q.empty()) {
        auto f = q.top();
        int edge = f.edge();
        q.pop();
        if (visit.count(edge) && visit[edge] >= f.it()) {
            continue;
        }

        // a terminal is found
        double w = 0.0;
        for (auto s: g.edge(edge).edge_signals()) {
            if (!signals.count(s)) {
                w += g.signal_weight(s);
            }
        }
        if (current[edge] && w > 0 && visit.find(edge) == visit.end()) {
            while(edge != -1 && in_solution.find(edge) == in_solution.end()) {
                q.emplace(edge, iteration, 0.0);
                shortest_distance[edge] = 0.0;
                in_solution.insert(edge);
                solution.add_edge(edge);
                for (auto s: g.edge(edge).edge_signals()) {
                    if (signals.find(s) == signals.end()) {
                        signals.insert(s);
                        solution.add_obj(g.signal_weight(s));
                    }
                }
                edge = previous[edge];
            }
            if (solution.objective() > best.objective()) {
                best = solution;
            }
        }

        edge = f.edge();
        iteration++;
        visit[edge] = iteration;

        // exploring all the neighbours
        int p = previous[edge];

        auto explore_neighbours = [this, edge, &q](size_t v) {
            for (auto& other: g.neighbours(v)) {
                if ((int)other.num() == edge) continue;
                double raise = std::max(1e-5, -weight(other.num()));
                double new_dist = shortest_distance[edge] + raise;
                if (new_dist < shortest_distance[other.num()]) {
                    previous[other.num()] = edge;
                    shortest_distance[other.num()] = new_dist;
                    q.emplace(other.num(), iteration, new_dist);
                }
            }
        };

        if (p == -1) {
            // first edge is considered
            explore_neighbours(g.edge(edge).from());
            explore_neighbours(g.edge(edge).to());
        } else {
            auto from = g.edge(edge).from();
            if (from == g.edge(p).from() || from == g.edge(p).to()) {
                explore_neighbours(g.edge(edge).to());
            } else {
                explore_neighbours(from);
            }
        }
    }
}

double PrimalHeuristic::wei(size_t e) {
    double local_wei = 0.0;
    for (auto s: g.edge(e).edge_signals()) {
        local_wei += g.signal_weight(s);
    }
    return local_wei;
}

}