//
// Created by Alexander Loboda on 9/10/21.
//

#ifndef MWCSR_COMPONENTS_H
#define MWCSR_COMPONENTS_H

#include "../../include/graph.h"

#include <vector>
#include <unordered_set>

namespace relax {

class Component {
    std::vector<size_t> edges;
    std::unordered_set<size_t> adjs;
    double max_revenue;
public:
    Component();
    Component(Component&& other) = default;
    Component& operator=(Component&& other) = default;
    void add_edge(size_t edge);
    void add_neighbour(size_t edge);
    bool operator<(const Component&) const;
    void set_revenue(double r);
    double get_revenue() const;

    std::vector<size_t> component_edges() const;
    std::vector<size_t> component_env() const;

    static std::vector<Component> get_components(mwcsr::Graph& g, std::vector<size_t> edges);
};

#endif //MWCSR_COMPONENTS_H

}
