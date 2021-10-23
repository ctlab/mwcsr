#ifndef ANNEALING_MODULE_H
#define ANNEALING_MODULE_H

#include "../../include/graph.h"

using namespace mwcsr;

namespace annealing {

    class Module {
        std::vector<size_t> vs;
        std::vector<Edge> es;

    public:
        Module();
        Module(const Graph& g, const std::vector<size_t>& vertices, const std::vector<size_t>& edges);

        std::vector<size_t> vertices();
        std::vector<Edge> edges();
    };

}

#endif //ANNEALING_MODULE_H
