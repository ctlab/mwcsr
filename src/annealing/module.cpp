#include "include/module.h"

namespace annealing {

    std::vector<size_t> Module::vertices() {
        return vs;
    }

    std::vector<Edge> Module::edges() {
        return es;
    }

    Module::Module() = default;

    Module::Module(const Graph& g, const std::vector<size_t>& vertices, const std::vector<size_t>& edges) {
        vs = vertices;
        for (size_t edge : edges) {
            es.push_back(g.const_edge(edge));
        }
    }
}
