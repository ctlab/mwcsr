//
// Created by Alexander Loboda on 9/10/21.
//

#include "include/components.h"

#include <unordered_set>
#include <algorithm>
#include <queue>

namespace relax {

Component::Component() {}

void Component::add_edge(size_t edge) {
    edges.push_back(edge);
}

void Component::add_neighbour(size_t edge) {
    adjs.insert(edge);
}

bool Component::operator<(const Component& other) const {
    return edges.size() < other.edges.size();
}

std::vector<Component> Component::get_components(mwcsr::Graph& g, std::vector<size_t> sol) {
    std::unordered_set<size_t> solution(sol.begin(), sol.end());
    std::unordered_set<size_t> visited;

    std::vector<Component> components;

    for (size_t start: sol) {
        if (visited.count(start)) {
            continue;
        }
        Component c;
        std::unordered_set<size_t> positive_signals;
        double max_revenue = 0.0;

        std::queue<size_t> q;
        visited.insert(start);
        q.push(start);
        while(!q.empty()) {
            auto e = q.front();
            q.pop();
            auto edge = g.edge(e);
            c.add_edge(edge.num());
            for (auto s: edge.edge_signals()) {
                if (g.signal_weight(s) > 0 && positive_signals.find(s) == positive_signals.end()) {
                    positive_signals.insert(s);
                    max_revenue += g.signal_weight(s);
                }
            }
            for (auto v: {g.edge(e).from(), g.edge(e).to()}) {
                for (auto e: g.neighbours(v)) {
                    if (!solution.count(e.num())) {
                        c.add_neighbour(e.num());
                        continue;
                    }

                    if (!(visited.find(e.num()) == visited.end())) {
                        continue;
                    }
                    visited.insert(e.num());

                    q.push(e.num());
                }
            }
        }

        c.set_revenue(max_revenue);
        components.push_back(std::move(c));
    }

    std::sort(components.begin(), components.end());
    return components;
}

std::vector<size_t> Component::component_edges() const {
    return edges;
}

std::vector<size_t> Component::component_env() const {
    return {adjs.begin(), adjs.end()};
}

void Component::set_revenue(double r) {
    max_revenue = r;
}

double Component::get_revenue() const {
    return max_revenue;
}

}
